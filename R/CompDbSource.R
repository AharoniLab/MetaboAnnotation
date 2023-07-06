#' @include CompAnnotationSource.R

#' @title Compound Annotation Sources for `CompDb` databases
#'
#' @aliases CompDbSource-class
#'
#' @description
#'
#' `CompDbSource` objects represent references to [CompDb] database-backed
#' annotation resources. Instances are expected to be created with the dedicated
#' construction functions such as `MassBankSource` or the generic
#' `CompDbSource`. The annotation data is not stored within the object but will
#' be accessed/loaded within the object's `matchSpectra` method.
#'
#' New `CompDbSource` objects can be created using the functions:
#'
#' - `CompDbSource`: create a new `CompDbSource` object from an existing
#'   `CompDb` database. The (SQLite) database file (including the full path)
#'   needs to be provided with parameter `dbfile`.
#'
#' - `MassBankSource`: retrieves a `CompDb` database for the specified MassBank
#'   release from Bioconductor's online `AnnotationHub` (if it exists) and
#'   uses that. Note that `AnnotationHub` resources are cached locally and thus
#'   only downloaded the first time.
#'   The function has parameters `release` which allows to define the desired
#'   MassBank release (e.g. `release = "2021.03"` or `release = "2022.06"`)
#'   and `...` which allows to pass optional parameters to the `AnnotationHub`
#'   constructor function, such as `localHub = TRUE` to use only the cached
#'   data and avoid updating/retrieving updates from the internet.
#'
#' Other functions:
#'
#' - `metadata`: get metadata (information) on the annotation resource.
#'
#' @param dbfile `character(1)` with the database file (including the full
#'     path).
#'
#' @param object A `CompDbSource` object.
#'
#' @param release A `character(1)` defining the version/release of MassBank that
#'     should be used.
#'
#' @param x A `CompDbSource` object.
#'
#' @param ... For `CompDbSource`: ignored. For `MassBankSource`: optional
#'     parameters passed to the `AnnotationHub` constructor function.
#'
#' @name CompDbSource
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Locate a CompDb SQLite database file. For this example we use the test
#' ## database from the `CompoundDb` package.
#' fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
#' ann_src <- CompDbSource(fl)
#'
#' ## The object contains only the reference/link to the annotation resource.
#' ann_src
#'
#' ## Retrieve a CompDb with MassBank data for a certain MassBank release
#' mb_src <- MassBankSource("2021.03")
#' mb_src
NULL

#' @description
#'
#' It is better to NOT put neither a connection object nor the data itself into
#' the source object to allow parallel processing. The `matchSpectra` method
#' needs then to make sure to load/retrieve the data and connect to it.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @exportClass CompDbSource
setClass(
    "CompDbSource",
    contains = "CompAnnotationSource",
    slots = c(dbfile = "character"),
    prototype = list(dbfile = character())
)

#' validator to check if the provided file name is indeed a CompDb database.
#'
#' @importFrom CompoundDb CompDb
#'
#' @noRd
.validate_dbfile <- function(x) {
    if (length(x))
        validObject(CompDb(x))
    else TRUE
}

setValidity("CompDbSource", function(object) {
    .validate_dbfile(object@dbfile)
})

#' @rdname CompDbSource
#'
#' @export
CompDbSource <- function(dbfile = character()) {
    new("CompDbSource", dbfile = dbfile)
}

#' @export
#'
#' @rdname CompDbSource
setMethod("metadata", "CompDbSource", function(x, ...) {
    db <- CompDb(x@dbfile)
    metadata(db)
})

#' @importFrom methods callNextMethod
#'
#' @rdname CompDbSource
setMethod("show", "CompDbSource", function(object) {
    callNextMethod()
    md <- metadata(object)
    cat("Metadata information:\n")
    cat(paste0("  - ", md[, 1], ": ", md[, 2], "\n"), sep = "")
})

#' @rdname CompareSpectraParam
#'
#' @importFrom Spectra MsBackendDataFrame
#'
#' @importClassesFrom Spectra MsBackendDataFrame
setMethod(
    "matchSpectra", signature(query = "Spectra", target = "CompDbSource",
                              param = "Param"),
    function(query, target, param, BPPARAM = BiocParallel::SerialParam()) {
        ## connect to the database
        db <- CompDb(target@dbfile)
        ## get the Spectra from the source and call matchSpectra
        res <- matchSpectra(query, Spectra(db), param = param,
                            BPPARAM = BPPARAM)
        ## keep only matching reference/target spectra  and change the
        ## backend to MsBackendDataFrame
        res <- pruneTarget(res)
        res@target <- setBackend(res@target, backend = MsBackendDataFrame())
        res
    })

#' @export
#'
#' @rdname CompDbSource
MassBankSource <- function(release = "2021.03", ...) {
    if (!requireNamespace("AnnotationHub", quietly = TRUE))
        stop("'MassBankSource' requires the 'AnnotationHub' package.\n",
             "Please install it with 'BiocManager::install(",
             "\"AnnotationHub\")' and try again.")
    ah <- AnnotationHub::AnnotationHub(...)
    res <- AnnotationHub::query(ah, c("MassBank", release))
    if (!length(res))
        stop("MassBank release \"", release, "\" not found in AnnotationHub")
    if (length(res) > 1)
        stop("Provided release is ambiguous: ", length(res), " data sets in ",
             "AnnotationHub match the provided release information.")
    fn <- unname(AnnotationHub::cache(res))
    CompDbSource(fn)
}

#' @title Annotation using MetFrag
#'
#' @rdname MetFrag
#'
#' @name MetFragParam,MetFragSource, validAdducts
#'
#' @description
#'
#' Experimental (query) fragment spectra can be matched against the MetFrag
#' server using `matchSpectra` by providing `MatchFragSource()` with parameter
#' `target`. The parameters for the matching against MetFrag can be specified
#' using a `MetFragParam` object (parameter `param` of `matchSpectra`).
#'
#' Description of functions:
#'
#' - `MetFragParam`: creates an `MetFragParam` instance with the parameters
#'   for the MetFrag query.
#'
#' - `MetFragSource`: creates a `MetFragSource` instance and checks if the
#'   specified URL is valid.
#'
#' - `matchSpectra`: annotates provided MS/MS spectra (as a [Spectra()] object
#'   with parameter `query`) against MetFrag. Polarity and expected adduct
#'   information have to be set in the provided `Spectra` object: polarity
#'   has to be encoded with `0` and `1` for negative and positive polarity
#'   through the spectra variable `"polarity"`. The expected adduct(s) for
#'   each fragment spectrum need to be defined as a spectra variable
#'   `"precursorAdduct"`. Use `validAdducts(MetFragParam())` to get the set of
#'   adduct definitions supported by MetFrag.
#'   Parameter `target` is expected to be an instance of `MetFragSource`
#'   providing the URL to the MetFrag server that should be used for the
#'   matching. Settings for the matching can be provided using an instance of
#'   `MetFragParam` with parameter `param`. The function iterates over all
#'   provided query spectra and annotates each using the MetFrag server. Note
#'   that an active internet connection is required for this function. Note
#'   also that queries can timeout (in which case a warning is shown and no
#'   match is reported for the specific fragment spectrum). The result is
#'   returned as a [Matched()] object with `query` being the input `Spectra`
#'   and `target` the `data.frame` with matches returned from the MetFrag
#'   server. See examples below for details.
#'
#' - `validAdducts`: `validAdducts(MetFragParam())` returns a `list` with the
#'   supported adduct names of MetFrag (for positive and negative polarity).
#'   If called on a [Spectra()] object (`validAdducts(spectra, MetFragParam())`)
#'   it checks whether the specified adduct names (with spectra variable
#'   `"precursorAdduct"`) are valid, given the polarity of each spectrum. The
#'   function returns `TRUE` if all adduct definitions are valid or throws an
#'   error.
#'
#' LLLLL PARAMETERS
#' @author Nir Shahaf, Johannes Rainer
NULL

## We Assume that the user has added the precursor adduct annotation(s) to the
## spectra data, under the column 'precursorAdduct' e.g.:
## load(system.file("extdata", "sps_test.Rdata",package = "MetaboAnnotation"))
## sps_test$precursorAdduct <- "[M-H]-"
## Define MetFrag parameters object and make connection with the online server:
## mfParam = MetFragParam()
## mfSrc = MetFragSource()
##
## testMatch <- matchSpectra(sps_test, MetFragSource(), MetFragParam())

#' @rdname MetFrag
#'
#' @export
setClass("MetFragParam",
         slots = c(
             validAdducts = "list",
             ## database parameters -> how to retrieve candidates
             metfragdatabasetype = "character",
             OfflineSpectralDatabaseFile = "character",
             ## peak matching parameters
             fragmentpeakmatchabsolutemassdeviation = "numeric",
             fragmentpeakmatchrelativemassdeviation = "numeric",
             databasesearchrelativemassdeviation = "numeric",
             ## scoring parameters
             MetFragScoreTypes = "character",
             MetFragScoreWeights = "numeric",
             ## default class filtering params
			 requirePrecursor = "logical",
             requirePrecursorPeak = "logical",
             THRESHFUN = "function",
             toleranceRt = "numeric",
             percentRt = "numeric",
			 minFragCount = "numeric"
         ),
         contains = "Param",
		 prototype = prototype(
			 validAdducts = list(
                 positive = c("[M]+","[M+H]+","[M+NH4]+","[M+Na]+","[M+K]+",
                              "[M+CH4O+H]+","[M+C2H3N+H]+","[M+C2H3N+Na]+",
                              "[M+C4H6N2+H]+"),
                 negative = c("[M]-","[M-H]-","[M+Cl]-","[M+CHO2]-",
                              "[M+C2H3O2]-")
			 ),
             metfragdatabasetype = "PubChem",
			 OfflineSpectralDatabaseFile =
                 "/vol/spectral-databases/weizfragV2.mb",
             fragmentpeakmatchabsolutemassdeviation = 0.0025,
             fragmentpeakmatchrelativemassdeviation = 15,
             databasesearchrelativemassdeviation = 7.5,
             MetFragScoreTypes = "FragmenterScore,AutomatedPeakFingerprintAnnotationScore",
             MetFragScoreWeights = c(1, 1),
			 requirePrecursor = TRUE,
             requirePrecursorPeak = FALSE,
             THRESHFUN = function(x) which(x >= 1),
             toleranceRt = Inf,
             percentRt = 0,
			 minFragCount = 3
         ),
         validity = function(object) {
             msg <- NULL
			 if (!all(unlist(object@validAdducts) %in%
                      c(adductNames("negative"), adductNames("positive"))
                      ))
                 msg <- paste0("'validAdducts' contains unrecognized adduct ",
                               "names. Check 'MetaboCoreUtils::adductNames' ",
                               "documentation")
             if (length(object@fragmentpeakmatchabsolutemassdeviation) > 1 ||
                 object@fragmentpeakmatchabsolutemassdeviation < 0)
                 msg <- c(msg,
                          paste0("'fragmentpeakmatchabsolutemassdeviation'",
                                 " has to be a positive number of length 1"))
             if (length(object@fragmentpeakmatchrelativemassdeviation) > 1 ||
                 object@fragmentpeakmatchrelativemassdeviation < 0)
                 msg <- c(msg,
                          paste0("'fragmentpeakmatchrelativemassdeviation' has",
                                 " to be positive number of length 1"))
              if (length(object@databasesearchrelativemassdeviation) > 1 ||
                  object@databasesearchrelativemassdeviation < 0)
                 msg <- c(msg,
                          paste0("'databasesearchrelativemassdeviation' has ",
                                 "to be positive number of length 1"))
             msg <- c(msg, .valid_threshfun(object@THRESHFUN))
             if (any(object@toleranceRt < 0))
                 msg <- c(msg, "'toleranceRt' has to be positive")
             if (any(object@percentRt < 0))
                 msg <- c(msg, "'percentRt' has to be positive")
			 if (object@minFragCount < 3)
				msg <- c(msg, paste0("Minimum number of input spectra ",
                                     "fragments needs to be at least three"))
             msg
         })

#' @rdname MetFrag
#'
#' @export
MetFragParam <- function(...) {
	new("MetFragParam")
}

#' @rdname MetFrag
#'
#' @export
setGeneric("validAdducts", function(object, ...)
    standardGeneric("validAdducts"))

#' @rdname MetFrag
#'
#' @export
setMethod("validAdducts", "MetFragParam", function (object,...)
    object@validAdducts)

#' @rdname MetFrag
setMethod("validAdducts", "Spectra",
          function(object, param = MetFragParam(), ...) {
              if (!any(spectraVariables(object) == "precursorAdduct"))
                  stop("'object' is expected to have a spectra variable ",
                       "\"precursorAdduct\" specifying which potential ",
                       "adduct(s) of the original molecule the precursor ",
                       "ion might be.")
              pol <- polarity(object)
              if (!all(pol %in% c(0, 1)))
                  stop("polarity information in 'object' is either missing or ",
                       "wrong. Polarity should be either 0 (negative) or 1 ",
                       "(positive).")
              ## Check positive
              if (!all(unlist(object$precursorAdduct[pol == 1]) %in%
                       validAdducts(param)$positive))
                  stop("Some of the specified adduct definitions (in spectra ",
                       "variable \"precursorAdduct\" are not valid. See ",
                       "`validAdducts(MetFragParam())` for supported adduct ",
                       "names ")
              if (!all(unlist(object$precursorAdduct[pol == 0]) %in%
                       validAdducts(param)$negative))
                  stop("Some of the specified adduct definitions (in spectra ",
                       "variable \"precursorAdduct\" are not valid. See ",
                       "`validAdducts(MetFragParam())` for supported adduct ",
                       "names ")
              TRUE
})

#' This class is defined by default to access the online server instance
#' running at ipb-halle: https://msbi.ipb-halle.de/MetFrag-deNBI/api/v1/process
#' the default server URL is set to the current address given by S.N. -
#' but this can be set to a local address running a docker container
#' Configuration options are as in the 'MetFragRelaunched_miniclient.R'
#' file suggested by S.N.
#'
#' @author Nir Shahaf
#'
#' @noRd
setClass(
    "MetFragSource",
    contains = "CompAnnotationSource",
    slots = c(url = "character",
              config = "list"),
    prototype = list(
        url = "https://msbi.ipb-halle.de/MetFrag-deNBI/api/v1/process",
        config = list(`Content-Type` = "application/json"))
)

#' @importFrom httr POST verbose
#'
#' @rdname MetFrag
#'
#' @export
MetFragSource <-
    function(url = "https://msbi.ipb-halle.de/MetFrag-deNBI/api/v1/process",
             debug = FALSE) {
        mfSrc <- new("MetFragSource", url = url)
        if (debug) handle <- verbose()
        else handle <- NULL
        testQuery <- list(
            fragmentpeakmatchabsolutemassdeviation = "0.001",
            fragmentpeakmatchrelativemassdeviation = "5",
            databasesearchrelativemassdeviation = "5",
            metfragdatabasetype = "PubChem",
            peakliststring = paste0("177.092320760091_1451.08410644531;",
                                    "205.086236000061_1983.02587890625;",
                                    "229.086613972982_1455.84973144531;",
                                    "234.125782775879_9424.7314453125;",
                                    "243.102464294434_6834.91943359375;",
                                    "244.109477996826_1533.8046875;",
                                    "244.146520514237_4892.56103515625;",
                                    "368.272043863932_1573.63305664062;",
                                    "369.206732177734_12602.3115234375;",
                                    "370.211002349854_1442.79187011719;",
                                    "371.222595214844_2184.82690429688;",
                                    "491.388732910156_4496.3466796875;",
                                    "492.324691772461_5253.52587890625"),
            neutralprecursormolecularformula = "C10H11N4O5",
            precursoriontype = "[M-H]-",
            neutralprecursormass = "267.07291",
            OfflineSpectralDatabaseFile =
                "/vol/spectral-databases/weizfragV2.mb"
        )
        testCall <- POST(url = mfSrc@url,config = mfSrc@config,
                         body = testQuery, encode = "json", handle)
        ## a naive way to capture http 'success' codes 2xx:
        if (length(testCall$status_code) &&
            length(grep("^2", testCall$status_code))) mfSrc
        else stop("Remote server inactive or irresponsive. ",
                  "Returned status code: ", testCall$status_code)
}

.peaks_to_metfrag_string <- function(x) {
    if (nrow(x))
        paste0(x[, "mz"], "_", x[, "intensity"], collapse = ";")
    else character()
}

#' @importMethodsFrom ProtGenerics as.list
.param_to_list <- function(x, skip = c("validAdducts", "THRESHFUN")) {
    x <- as.list(x)
    x[!names(x) %in% skip]
}

#' @rdname MetFrag
#'
#' @importFrom httr GET POST content
#'
#' @importFrom utils read.csv
#'
#' @importFrom MetaboCoreUtils mz2mass
setMethod(
	"matchSpectra", signature(query = "Spectra", target = "MetFragSource",
                              param = "MetFragParam"),
	function(query, target, param, debug = FALSE) {
        validAdducts(query, param)
        if (debug) handle <- verbose()
        else handle <- NULL
        out <- vector("list", length(query))
		names(out) <- query$peak_id
        restQuery <- .param_to_list(param)
		for (i in seq_along(query)) {
            qi <- query[i]
			peaks <- peaksData(qi)[[1L]]
			if (nrow(peaks) < param@minFragCount) {
				## In cases of very few fragment peaks, the MetFrag server is
                ## 'timing-out' due (probably) to the large search space -
				##	thus, I think the it's better to send the query with an
                ## empty fragment spectra and at least get the basic initial
                ## search matrix.
				## ...the other option is to ignor these calls - leaving the
                ## input as-is and returning empty results:
				# out[[i]] = list(
					# mfResults = NULL,
					# mfServerStatus = "NA"
				# )
				# next()
                peakListString <- ""
			} else
                peakListString <- .peaks_to_metfrag_string(peaks)
			pmz <- qi$precursorMz
			adduct <- unlist(qi$precursorAdduct)
			out[[i]] <- do.call(rbind, lapply(adduct, function(a) {
				query <- c(
					restQuery,
					neutralprecursormass = as.character(mz2mass(pmz, a)),
					precursoriontype = a,
					peakliststring = peakListString
				)
				mfRequest <- POST(url = target@url, config = target@config,
                                  body = query, encode = "json", handle)
				## get the http call status, plus results(if exist):
				callStatus <- .chk_rest_response(mfRequest)
				if (callStatus$status == "SUCCESS" ) {
					res <- read.csv(text = content(GET(callStatus$rest_url)))
                    cbind(res, data.frame(query_idx = rep(i, nrow(res)),
                                          adduct = rep(a, nrow(res))))
                } else {
                    warning("No result for query spectrum [", i, "], adduct \"",
                            a, "\". Response from server was: ",
                            callStatus$status)
                    ## Suggestion: don't report anything if MetFrag did not
                    ## return any result (no matter why). Alternative would be
                    ## to report `NA` for all fields, which might however be
                    ## misleading to the user (because a *match* would be
                    ## reported).
                    data.frame()
                }
			}))
			Sys.sleep(1)
		}
		out <- do.call(rbind, out)
        mtches <- data.frame(query_idx = out$query_idx,
                             target_idx = seq_len(nrow(out)),
                             score = out$Score,
                             adduct = out$adduct)
        Matched(query = query,
                target = out[, !colnames(out) %in%
                               c("query_idx", "Score", "adduct")],
                matches = mtches, metadata = list(param))
	}
)

.chk_rest_response <- function(response, gracetime = 2, max_try = 10, ...) {
	statusIdx <- which(names(content(response)$"_links") == "status")
	STATUSURL <- content(response)$"_links"[[statusIdx]]$href
	STATUSURL <- sub ('^http','https', STATUSURL)[[1]]

	resultIdx <- which(names(content(response)$"_links") == "result")
	RESULTURL <- content(response)$"_links"[[resultIdx]]$href
	RESULTURL <- sub ('^http','https',RESULTURL)[[1]]

	for (i in seq_len(max_try)) {
		Sys.sleep(gracetime)
		gracetime <- gracetime + i
		status <- content(GET(STATUSURL))$status
		if (status == "RUNNING") {
			cat(".")
			status  <- "TIMEOUT" ## If this is not the last iteration, TIMEOUT will be overwritten by next status poll
			next
		} else if (status == "SUCCESS") {
			cat("Status: ", status, "\n")
			break
		} else {
            break
		}
	}
	list(status = status, rest_url = RESULTURL)
}
