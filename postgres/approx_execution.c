/**********************************************************************
 *
 * FESTIval - Framework to Evaluate SpaTial Indices in non-VolAtiLe memories and hard disk drives.
 * https://accarniel.github.io/FESTIval/
 *
 * Copyright (C) 2016-2020 Anderson Chaves Carniel <accarniel@gmail.com>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU General Public Licence. See the COPYING file.
 *
 * Fully developed by Anderson Chaves Carniel
 *
 **********************************************************************/

/* this is the main file that define the postgresql functions for spatial approximations */
/*needed libraries from postgres*/
#include "postgres.h"
#include "funcapi.h"
#include "fmgr.h"
#include "executor/spi.h"
#include "utils/memutils.h"

#include <math.h>
#include <liblwgeom.h>
#include <unistd.h>    //to handle POSTGIS objects
#include "lwgeom_pg.h" //to convert POSTGIS objects

#include "../main/spatial_index.h"
#include "../main/log_messages.h"
#include "../main/statistical_processing.h"
#include "../main/lwgeom_handler.h"
#include "../approximations/SpatialApproximation.h"
#include "../approximations/SpatialApproximation_utils.h"

#include "query.h" //for query processing

#include "access/htup_details.h"
#include "utils/builtins.h"

#include "../festival_config.h"

static Source *read_source_from_fds(int src_id);

Source *read_source_from_fds(int src_id)
{
    char query[256];
    int err;
    MemoryContext old_context;

    char *p, *t, *c, *s;

    Source *src = (Source *)lwalloc(sizeof(Source));

    sprintf(query, "SELECT schema_name, table_name, column_name, pk_name "
                   "FROM fds.source WHERE src_id = %d;",
            src_id);

    if (SPI_OK_CONNECT != SPI_connect())
    {
        SPI_finish();
        _DEBUG(ERROR, "read_source_from_fds: could not connect to SPI manager");
        return NULL;
    }
    err = SPI_execute(query, true, 1);
    if (err < 0)
    {
        SPI_finish();
        _DEBUG(ERROR, "read_source_from_fds: could not execute the SELECT command");
        return NULL;
    }

    if (SPI_processed <= 0)
    {
        SPI_finish();
        _DEBUGF(ERROR, "the src_id (%d) does not exist in the table", src_id);
        return NULL;
    }

    s = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 1);
    t = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 2);
    c = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 3);
    p = SPI_getvalue(SPI_tuptable->vals[0], SPI_tuptable->tupdesc, 4);

    /* get the variables in upper memory context (outside of SPI) */
    /* TODO: use a narrower context to switch to */
    old_context = MemoryContextSwitchTo(TopMemoryContext);
    src->schema = (char *)lwalloc(strlen(s) + 1);
    strcpy(src->schema, s);

    src->table = (char *)lwalloc(strlen(t) + 1);
    strcpy(src->table, t);

    src->column = (char *)lwalloc(strlen(c) + 1);
    strcpy(src->column, c);

    src->pk = (char *)lwalloc(strlen(p) + 1);
    strcpy(src->pk, p);

    src->src_id = src_id;

    MemoryContextSwitchTo(old_context);

    pfree(s);
    pfree(t);
    pfree(c);
    pfree(p);

    /* disconnect from SPI */
    SPI_finish();

    return src;
}

/*  */
Datum STI_build_approximations(PG_FUNCTION_ARGS);

PG_FUNCTION_INFO_V1(STI_build_approximations);

Datum STI_build_approximations(PG_FUNCTION_ARGS)
{
    //args and pg variables
    text *path = PG_GETARG_TEXT_P(0);
    int approx_id = PG_GETARG_INT32(1);
    int src_id = PG_GETARG_INT32(2);
    Source *src;
    char *file_path;

    //query essentials
    char query[256];
    int err;
    bytea *bytea_ewkb;
    uint8_t *ewkb;
    char isnull;
    char *ptemp;

    //spatial approx. builder variables
    size_t approx_size;
    LWGEOM **geoms;
    int *gids;
    SpatialApproximation **approx;
    const SpatialApproximation **ap;
    FileSpecification fs;
    int num_per_page, num_per_query, total, cur_offset, i, j, k, c_max, increment, acum;

    MemoryContext oldercontext, oldcontext;
    // this will survives until a restart of the server
    oldercontext = MemoryContextSwitchTo(TopMemoryContext);

#if FESTIVAL_PGSQL_VERSION >= 120

    file_path = text_to_cstring(path);

#else

    file_path = text2cstring(path);

#endif

    if (approx_id != BBOX_TYPE && approx_id != MBC_TYPE && approx_id != RMBR_TYPE && approx_id != RMBP_TYPE && approx_id != MBE_TYPE && approx_id != N_CORNER_5_TYPE && approx_id != N_CORNER_4_TYPE && approx_id != 7)
    {
        ereport(ERROR,
                (errcode(ERRCODE_CASE_NOT_FOUND),
                 errmsg("There is no approximation associated with the id provided - %d", approx_id)));
    }

    if (approx_id == 7)
    {
        approx_id = N_CORNER_5_TYPE;
    }

    src = read_source_from_fds(src_id);

    /* setting variables */
    approx_size = spatialapproximation_size(approx_id);

    fs.index_path = lwalloc(strlen(file_path) + strlen(get_approx_name(approx_id)) + 2);
    sprintf(fs.index_path, "%s.%s", file_path, get_approx_name(approx_id));
    fs.page_size = 4096;
    fs.io_access = DIRECT_ACCESS;

    num_per_page = floor(((double)fs.page_size - sizeof(uint32_t)) / (double)approx_size);
    num_per_query = num_per_page * 20;
    cur_offset = 0;
    c_max = 0;
    increment = 0;
    acum = 0;
    ap = lwalloc(num_per_page * sizeof(SpatialApproximation *));

    /* Query and file I/O */
    do
    {
        sprintf(query, "SELECT %s, ST_AsEWKB(%s) FROM %s.%s ORDER BY %s LIMIT %d OFFSET %d;", src->pk, src->column, src->schema, src->table, src->pk, num_per_query, cur_offset);

        if (SPI_OK_CONNECT != SPI_connect())
        {
            SPI_finish();
            _DEBUG(ERROR, "STI_build_approximations: could not connect to SPI manager");
            PG_RETURN_BOOL(false);
        }
        err = SPI_execute(query, true, 0);
        if (err < 0)
        {
            SPI_finish();
            _DEBUG(ERROR, "STI_build_approximations: could not execute the EXECUTE command");
            PG_RETURN_BOOL(false);
        }

        if (SPI_processed <= 0)
        {
            SPI_finish();
            _DEBUG(ERROR, "STI_build_approximations: returned 0 tuples");
            PG_RETURN_BOOL(false);
        }

        total = SPI_processed;

        oldcontext = MemoryContextSwitchTo(TopMemoryContext);

        geoms = (LWGEOM **)lwalloc(sizeof(LWGEOM *) * total);
        gids = lwalloc(total * sizeof(int));

        for (i = 0; i < total; i++)
        {

#if FESTIVAL_PGSQL_VERSION == 95

            bytea_ewkb = DatumGetByteaP(SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 2, &isnull));

#elif FESTIVAL_PGSQL_VERSION >= 120

            bytea_ewkb = DatumGetByteaP(SPI_getbinval(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 2, (bool *)&isnull));

#endif
            ewkb = (uint8_t *)VARDATA(bytea_ewkb);
            geoms[i] = lwgeom_from_wkb(ewkb, VARSIZE(bytea_ewkb) - VARHDRSZ, LW_PARSER_CHECK_NONE);

            ptemp = SPI_getvalue(SPI_tuptable->vals[i], SPI_tuptable->tupdesc, 1);
            gids[i] = atoi(ptemp);
            pfree(ptemp);
        }

        MemoryContextSwitchTo(oldcontext);
        /* disconnect from SPI */
        SPI_finish();

        approx = spatialapproximation_builder(total, geoms, approx_id, gids);
        acum = 0;
        for (i = 0; i < total; i += num_per_page)
        {
            acum += num_per_page;
            c_max += num_per_page;

            if (acum > total)
            {
                c_max -= (acum - total);
            }

            increment = c_max - i;

            k = 0;

            for (j = i; j < c_max; j++)
            {
                ap[k] = approx[j];
                k++;
            }

            spatialapproximation_write_to_file(approx_id, &fs, ap, increment);
        }

        cur_offset += num_per_query;
        c_max = 0;

        spatialapproximation_free_array(total, approx);
        lwgeom_free_array(total, geoms);
        lwfree(gids);
    } while (total == num_per_query);

    lwfree(ap);
    lwfree(fs.index_path);

    MemoryContextSwitchTo(oldercontext);
    PG_RETURN_BOOL(true);
}

/*index_name, index_path, type_query, query_bbox, predicate, approx_paths[]*/
PG_FUNCTION_INFO_V1(STI_query_spatial_index_with_approx);

Datum STI_query_spatial_index_with_approx(PG_FUNCTION_ARGS)
{

#if FESTIVAL_PGSQL_VERSION >= 120

    char *index_name = text_to_cstring(PG_GETARG_TEXT_PP(0));
    char *index_path = text_to_cstring(PG_GETARG_TEXT_PP(1));

#else

    char *index_name = text2cstring(PG_GETARG_TEXT_PP(0));
    char *index_path = text2cstring(PG_GETARG_TEXT_PP(1));

#endif

    int type_query = PG_GETARG_INT32(2);
    LWGEOM *lwgeom;
    GSERIALIZED *geom = PG_GETARG_GSERIALIZED_P(3);
    int predicate = PG_GETARG_INT32(4);
    int type_of_processing = FILTER_REFINEMENT_AND_APPROX_STEPS;
    char *spc_path;

    SpatialIndex *si;
    QueryResult *result;

    ReturnSetInfo *rsinfo = (ReturnSetInfo *)fcinfo->resultinfo;
    Tuplestorestate *tupstore;
    TupleDesc tupdesc;
    uint64 call_cntr;
    uint64 max_calls;
    MemoryContext per_query_ctx;
    MemoryContext oldcontext;

    ArrayType *array;
    int nelems;
    FileSpecification *approx_paths;
    int nof_approx;

    ArrayIterator iterator;
    Datum value;
    bool isnull;

    array = PG_GETARG_ARRAYTYPE_P(5);
    nelems = ArrayGetNItems(ARR_NDIM(array), ARR_DIMS(array));
    
    /* check to see if caller supports us returning a tuplestore */
    if (rsinfo == NULL || !IsA(rsinfo, ReturnSetInfo))
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("set-valued function called in context that cannot accept a set")));
    if (!(rsinfo->allowedModes & SFRM_Materialize))
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("materialize mode required, but it is not "
                        "allowed in this context")));

    per_query_ctx = rsinfo->econtext->ecxt_per_query_memory;

    // this will survives until a restart of the server in order to collect statistical data
    oldcontext = MemoryContextSwitchTo(TopMemoryContext);
    lwgeom = lwgeom_from_gserialized(geom);

    /*checking if the input geometry is valid*/
    if (lwgeom_is_empty(lwgeom))
    {
        _DEBUG(ERROR, "This is an empty geometry");
    }
    if (type_query == POINT_QUERY_TYPE)
    {
        if (lwgeom->type != POINTTYPE)
        {
            _DEBUGF(ERROR, "Invalid geometry type (%d) for the POINT_QUERY_TYPE", lwgeom->type);
        }
    }
    else if (type_query == RANGE_QUERY_TYPE)
    {
        //we have to check if it is a rectangle, otherwise we have to convert it to its BBOX
        if (lwgeom->type == POLYGONTYPE)
        {
            //we check if it is rectangular-shaped
            LWPOLY *poly = lwgeom_as_lwpoly(lwgeom);
            //here, we only check the number of rings (which must be equal to 1) and its number of points (=5)
            if (poly->nrings != 1 && poly->rings[0]->npoints != 5)
            {
                _DEBUG(ERROR, "Invalid geometry format for RANGE_QUERY_TYPE");
                //TO-DO perhaps we can have a better checker
            }
        }
        else
        {
            //we consider its bbox
            LWGEOM *input;
            BBox *bbox = bbox_create();
            if ((!lwgeom->bbox))
            {
                lwgeom_add_bbox(lwgeom);
            }
            gbox_to_bbox(lwgeom->bbox, bbox);
            input = bbox_to_geom(bbox);
            lwgeom_free(lwgeom);
            lwgeom = input;
            lwfree(bbox);
        }
    }

    spc_path = lwalloc(strlen(index_name) + strlen(index_path) + strlen(".header") + 1);
    strcpy(spc_path, index_path);
    strcat(spc_path, index_name);
    strcat(spc_path, ".header");

    si = spatialindex_from_header(spc_path);

    /* setting up FileSpecifications from Array */

    approx_paths = lwalloc(sizeof(FileSpecification) * nelems);
    nof_approx = 0;
    iterator = array_create_iterator(array, 0, NULL);

    while (array_iterate(iterator, &value, &isnull))
    {
        text *tmp_text;
        char *tmp_char;

        FileSpecification fs;

        if (isnull)
            continue;
            
        tmp_text = (text *)DatumGetTextPP(value);
        
#if FESTIVAL_PGSQL_VERSION >= 120

        tmp_char = text_to_cstring(tmp_text);
        
#else

        tmp_char = text2cstring(tmp_text);

#endif

        fs.io_access = DIRECT_ACCESS;
        fs.page_size = 4096;
        fs.index_path = lwalloc(strlen(tmp_char)+1);
        strcpy(fs.index_path, tmp_char);
        
        lwfree(tmp_char);

        approx_paths[nof_approx] = fs;
        nof_approx++;
    }

    /*the index_time is collected inside this function
     the reason is that the query is processed in two steps: filtering and refinement*/

    result = process_spatial_selection(si, lwgeom, predicate, type_query, type_of_processing, approx_paths, nof_approx);

    lwgeom_free(lwgeom);
    lwfree(spc_path);
    lwfree(index_name);
    lwfree(index_path);
    lwfree(approx_paths);

    //back to the current context
    MemoryContextSwitchTo(oldcontext);

    /*connect to SPI manager*/
    if (SPI_OK_CONNECT != SPI_connect())
    {
        _DEBUG(ERROR, "could not connect to SPI manager");
    }

    /* get a tuple descriptor for our result type */
    switch (get_call_result_type(fcinfo, NULL, &tupdesc))
    {
    case TYPEFUNC_COMPOSITE:
        /* success */
        break;
    case TYPEFUNC_RECORD:
        /* failed to determine actual type of RECORD */
        ereport(ERROR,
                (errcode(ERRCODE_FEATURE_NOT_SUPPORTED),
                 errmsg("function returning record called in context "
                        "that cannot accept type record")));
        break;
    default:
        /* result type isn't composite */
        ereport(ERROR,
                (errcode(ERRCODE_DATATYPE_MISMATCH),
                 errmsg("return type must be a row type")));
        break;
    }

    /* switch to long-lived memory context
     */
    oldcontext = MemoryContextSwitchTo(per_query_ctx);

    /* make sure we have a persistent copy of the result tupdesc */
    tupdesc = CreateTupleDescCopy(tupdesc);

    /* initialize our tuplestore in long-lived context */
    tupstore = tuplestore_begin_heap(rsinfo->allowedModes & SFRM_Materialize_Random,
                                     false, 1024);

    MemoryContextSwitchTo(oldcontext);

    /* total number of tuples to be returned */
    max_calls = result->nofentries;

    for (call_cntr = 0; call_cntr < max_calls; call_cntr++)
    {
        HeapTuple tuple;

        Datum values[2];
        bool *nulls;

        /* build the tuple and store it */
        nulls = palloc(tupdesc->natts * sizeof(bool));
        nulls[0] = false;

        values[0] = Int32GetDatum(result->row_id[call_cntr]);
        if ((type_of_processing == FILTER_AND_REFINEMENT_STEPS || type_of_processing == FILTER_REFINEMENT_AND_APPROX_STEPS) && result->geoms[call_cntr])
        {
            /*setting this column to be NOT NULL(It is very important)*/
            values[1] = PointerGetDatum(geometry_serialize(result->geoms[call_cntr]));
            nulls[1] = false;
        }
        else
        {
            values[1] = (Datum)0;
            nulls[1] = true;
        }

        tuple = heap_form_tuple(tupdesc, values, nulls);
        pfree(nulls);
        tuplestore_puttuple(tupstore, tuple);

        heap_freetuple(tuple);
    }

    /* let the caller know we're sending back a tuplestore */
    rsinfo->returnMode = SFRM_Materialize;
    rsinfo->setResult = tupstore;
    rsinfo->setDesc = tupdesc;

    query_result_free(result, type_of_processing);

    PG_FREE_IF_COPY(geom, 3);
    /* release SPI related resources (and return to caller's context) */
    SPI_finish();

    return (Datum)0;
}