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

#' @name MetFragParam,MetFragSource, validAdducts
#'
#' @author Nir Shahaf, Johannes Rainer
#'
#' @examples

#' We Assume that the user has added the precursor adduct annotation(s) to the spectra data,
#' 	under the column 'precursorAdduct' e.g.:
#' load(system.file("extdata", "sps_test.Rdata",package = "MetaboAnnotation"))
#' spectraData(sps_test)$precursorAdduct = rep("[M-H]-", nrow = length(sps_test))
#' Define MetFrag parameters object and make connection with the online server:
#' mfParam = MetFragParam()
#' mfSrc = MetFragSource()
#' 
#' testMatch <- matchSpectra(sps_test,mfSrc,MetFragParam())
setClass("MetFragParam",
         slots = c(
		 #  validAdducts
			validAdducts = "list",
		 # database parameters -> how to retrieve candidates
			metfragdatabasetype = "character",
			OfflineSpectralDatabaseFile = "character",
		 # peak matching parameters	
			fragmentpeakmatchabsolutemassdeviation = "numeric",
			fragmentpeakmatchrelativemassdeviation = "numeric",
			databasesearchrelativemassdeviation = "numeric",
		 # scoring parameters
			MetFragScoreTypes = "character",
			MetFragScoreWeights = "numeric",
		 # default class filtering params 	
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
				"positive" = c("[M]+","[M+H]+","[M+NH4]+","[M+Na]+","[M+K]+","[M+CH4O+H]+","[M+C2H3N+H]+","[M+C2H3N+Na]+","[M+C4H6N2+H]+"),
				"negative" = c("[M]-","[M-H]-","[M+Cl]-","[M+CHO2]-","[M+C2H3O2]-")
			 ),
             metfragdatabasetype = "PubChem",
			 OfflineSpectralDatabaseFile = "/vol/spectral-databases/weizfragV2.mb",
             fragmentpeakmatchabsolutemassdeviation = 0.0025,
             fragmentpeakmatchrelativemassdeviation = 15,
             databasesearchrelativemassdeviation = 7.5,
             MetFragScoreTypes = "FragmenterScore,AutomatedPeakFingerprintAnnotationScore",
             MetFragScoreWeights = c(1,1),
			 requirePrecursor = TRUE,
             requirePrecursorPeak = FALSE,
             THRESHFUN = function(x) which(x >= 1),
             toleranceRt = Inf,
             percentRt = 0,
			 minFragCount = 3
         ),
         validity = function(object) {
             msg <- NULL
			 if (!all (
					unlist(object@validAdducts) %in% c(MetaboCoreUtils::adductNames("negative"),MetaboCoreUtils::adductNames("positive"))
				))
				msg <- c("'validAdducts' contains unrecognized adduct names. Check 'MetaboCoreUtils::adductNames' documentation")
             if (length(object@fragmentpeakmatchabsolutemassdeviation) > 1 || object@fragmentpeakmatchabsolutemassdeviation < 0)
                 msg <- c("'fragmentpeakmatchabsolutemassdeviation' has to be a positive number of length 1")
             if (length(object@fragmentpeakmatchrelativemassdeviation) > 1 || object@fragmentpeakmatchrelativemassdeviation < 0)
                 msg <- c("'fragmentpeakmatchrelativemassdeviation' has to be positive number of length 1")
              if (length(object@databasesearchrelativemassdeviation) > 1 || object@databasesearchrelativemassdeviation < 0)
                 msg <- c("'databasesearchrelativemassdeviation' has to be positive number of length 1")
             msg <- c(msg, .valid_threshfun(object@THRESHFUN))
             if (any(object@toleranceRt < 0))
                 msg <- c("'toleranceRt' has to be positive")
             if (any(object@percentRt < 0))
                 msg <- c("'percentRt' has to be positive")
			 if (object@minFragCount < 3)
				msg <- c("Minimum number of input spectra fragments needs to be at least three")
             msg
         })

#' Default constructor function for 'MetFragParam'
MetFragParam <- function(...) {

	new("MetFragParam")
	
}

#' Define validation calls for the precursor peak adduct annotations:
setGeneric("validAdducts", function(object, ...) standardGeneric("validAdducts"))

setMethod("validAdducts", "MetFragParam", function (object,...) object@validAdducts)

setMethod("validAdducts", "Spectra", function (object,param = MetFragParam(),...) 
{ 
	#' check polarity from spectra
	X = spectraData(object)	
	stopifnot("precursorAdduct" %in% names(X))
	stopifnot(all(X$polarity %in% c(0,1)))
	#' Match to corresponding adduct ector 
	validateAdducts = apply(X,1,function(x) {
	
		polarity = ifelse(x["polarity"] == 0,"negative","positive")
		x["precursorAdduct"] %in% validAdducts(MetFragParam())[[polarity]]
		
	})
	#' return or stop
	if (!all(validateAdducts)) {
		warning("Not all precursor adducts are valid! See details above")
		print(X[which(validateAdducts == FALSE),c("precursorMz","precursorAdduct")])
	} 
	# object@validAdducts)
	return(validateAdducts)
})

#' 'MetFragSource'
#' This class is defined by default to access the online server instance running at ipb-halle: '"https://msbi.ipb-halle.de/MetFrag-deNBI/api/v1/process"'
#' 	the default server URL is set to the current address given by S.N. - but this can be set to a local address running a docker container 
#' 	Configuration options are as in the 'MetFragRelaunched_miniclient.R' file suggested by S.N.  
setClass(
    "MetFragSource",
    contains = "CompAnnotationSource",
    slots = c(url = "character",
              config = "list"),    
    prototype = list(url = "https://msbi.ipb-halle.de/MetFrag-deNBI/api/v1/process", 
                     config = list(`Content-Type` = "application/json"))
)

#' Define the default constructor function for 'MetFragSource' and run a call to the server using a test query
MetFragSource <- function(url = "https://msbi.ipb-halle.de/MetFrag-deNBI/api/v1/process") {
    mfSrc <- new("MetFragSource", url = url)
	testQuery <- list(
		"fragmentpeakmatchabsolutemassdeviation" = "0.001",
		"fragmentpeakmatchrelativemassdeviation" = "5",
		"databasesearchrelativemassdeviation" = "5",
		"metfragdatabasetype" = "PubChem",
		"peakliststring" = "177.092320760091_1451.08410644531;205.086236000061_1983.02587890625;229.086613972982_1455.84973144531;234.125782775879_9424.7314453125;243.102464294434_6834.91943359375;244.109477996826_1533.8046875;244.146520514237_4892.56103515625;257.117561340332_1164.64953613281;261.112461853027_3875.52514648438;262.119772774833_1165.22351074219;271.133573913574_13089.0126953125;272.139910529642_2833.46997070312;274.120386297053_1286.1396484375;275.12854309082_28140.5546875;276.132524490356_2414.30102539062;285.14973449707_1103.77197265625;287.128471374512_2708.40209960938;287.164703369141_1599.52160644531;299.165280151367_1428.58227539062;300.136039733887_1403.30676269531;311.165107727051_1970.57592773438;313.149806213379_9827.076171875;313.180361938477_197179.75;314.152227228338_1240.31591796875;314.183631896973_14364.0771484375;315.159642028809_287124;316.163235473633_21380.60546875;325.180702209473_1367.2001953125;327.196063995361_1536.09033203125;329.175514622738_2210.23046875;330.182575426604_2943.61254882812;339.197692871094_1151.18811035156;340.241302490234_1050.67346191406;341.175484793527_1900.67651367188;343.190165201823_1245.60681152344;353.210777282715_1409.009765625;353.248142496745_3203.60986328125;354.25524597168_5494.884765625;355.191353352865_1707.62768554688;355.263327026367_6431.8408203125;367.263783772786_1120.03723144531;368.272043863932_1573.63305664062;369.206732177734_12602.3115234375;370.211002349854_1442.79187011719;371.222595214844_2184.82690429688;381.205629047595_6095.37890625;381.243681335449_13951.033203125;382.248254949396_1968.34484863281;383.222082519531_340943.9375;384.228602600098_73401.9609375;385.234485202365_3912.12866210938;395.223587036133_1795.80224609375;395.258506774902_73514.359375;396.263960266113_10231.8681640625;397.23766784668_204010.359375;398.242947387695_23262.927734375;399.249359130859_1100.38000488281;409.23782602946_2197.42846679688;410.245188395182_3128.14868164062;411.253760086863_8162.02685546875;422.31834763747_1995.22106933594;423.252960205078_24539.48046875;423.297088623047_1005.40423583984;424.257669448853_3149.98852539062;438.313815030185_1848.50830078125;451.284339904785_40680.7421875;452.28828671104_3518.17065429688;464.32948811849_1396.36779785156;465.29237874349_6912.3564453125;465.34294128418_2235.52661132812;466.307887268066_1091252.5;467.311297607422_100971.46875;491.388732910156_4496.3466796875;492.324691772461_5253.52587890625",
		"neutralprecursormolecularformula" = "C10H11N4O5",
		"precursoriontype" = "[M-H]-", 
		"neutralprecursormass" = "267.07291",
		"OfflineSpectralDatabaseFile" = "/vol/spectral-databases/weizfragV2.mb"	
	)
	testCall <- POST(url=mfSrc@url,config = mfSrc@config,body = testQuery,encode = "json",verbose())
	## a naive way to capture http 'success' codes 2xx:
	if (grep("^2",testCall$status_code) == 1) { return(mfSrc) }
	else { warning(paste ("Remote sever inactive or else irresponsive and returned the following statuse code:"),testCall$status_code); return(NULL) }
}

setMethod(
	"matchSpectra", signature(query = "Spectra", target = "MetFragSource",
                              param = "MetFragParam"),
	function(query, target, param, BPPARAM = BiocParallel::SerialParam()) {
        out <- list()
		names(out) <- spectraData(query)$peak_id
        mtch <- data.frame()	
		for (i in seq_along(query)) {
			## get peaks data from query[i] as list/character string for MetFrag (assuming single peaks table per spectra - is it always valid?!)	
			peaks <- peaksData(query[i])[[1]]
			## We allow either 'minFragCount' number of fragment peaks in spectra (under important uncontrolled auumption of proper data pretreatment), OR 
			##	a 
			if (nrow(peaks) < mfParam@minFragCount & nrow(peaks) > 0) {
				## In cases of very few fragment peaks, the MetFrag server is 'timing-out' due (probably) to the large search space -
				##	thus, I think the it's better to send the query with an empty fragment spectra and at least get the basic initial search matrix.
				peaks <- peaks[-c(1:nrow(peaks)),]
				## ...the other option is to ignor these calls - leaving the input as-is and returning empty results:
				# out[[i]] = list(
					# mfResults = NULL,
					# mfServerStatus = "NA"
				# )
				# next()				
			}
			
			peakListString <- paste(apply(peaks, 1, function (x){ paste(x[1],x[2],sep = '_') }),collapse = ';')

			mz <- spectraData(query[i])$precursorMz

			chkAdducts <- validAdducts(query,param)	
			
			if (!all(chkAdducts)) { 
				stop("Some spectra precursor adducts are not legal MetFrag adduct types: Use function 'validAdducts(MetFragParam())' to print a list of valid adduct names. 'matchSpectra' stopping.") 
			}
							
			adduct <- spectraData(query)$precursorAdduct[i]
			## All REST calls will be done inside this loop here:
			out[[i]] <- lapply(1:length(adduct),function (i) {
				neutralPrecursorMass <- MetaboCoreUtils::mz2mass(mz, adduct[i])
				restQuery <- sapply(slotNames(param)[-which(slotNames(param) %in% c("validAdducts","THRESHFUN"))], function(x) { 
					slot(param,x)
				})
				restQuery <- c(
					restQuery,	
					"neutralprecursormass" = as.character(neutralPrecursorMass),
					"precursoriontype" = as.character(adduct),
					"peakliststring" = peakListString
				)
				
				mfRequest <- POST(url = target@url,config = target@config,body = restQuery,encode = "json",verbose())
						
				## get the http call status, plus results(if exist): 
				callStatus <- .chk_rest_response (mfRequest)
				if (callStatus$status == "SUCCESS" ) {
					mfCandidates <- read.csv(text = content(GET(callStatus$rest_url))) 
				} else {				
					mfCandidates = NULL
					## Warning msg?
				}
				## return a list of anntation results table and server staus code
				list(
					mfResults = mfCandidates,
					mfServerStatus = callStatus$status
				)
			})
			names(out[[i]]) <- adduct
			Sys.sleep(1)
		}
		out	
		#' ToDo:
		#' Convert nested list to the MatchedSpectra Object?
		#' ...
		#' ...		
	}
)

.chk_rest_response <- function(response,gracetime = 2,...) {

	statusIdx <- which(names(content(response)$"_links") == "status")
	STATUSURL <- content(response)$"_links"[[statusIdx]]$href
	STATUSURL <- sub ('^http','https',STATUSURL)[[1]]

	resultIdx <- which(names(content(response)$"_links") == "result")
	RESULTURL <- content(response)$"_links"[[resultIdx]]$href
	RESULTURL <- sub ('^http','https',RESULTURL)[[1]]

	for (i in 1:10) {    
		
		Sys.sleep(gracetime)
		gracetime <- gracetime+i

		status <- content(GET(STATUSURL))$status

		if ( status == "RUNNING" ) {
			cat(".")
			status  <- "TIMEOUT" ## If this is not the last iteration, TIMEOUT will be overwritten by next status poll
			next
		} else if ( status == "SUCCESS" ) {
			cat("Status: ", status, "\n")
			break
		} else {
			cat("Status: ", status)
			cat("Error on server side...")
		}
	}
	return(list (status = status, rest_url = RESULTURL))
}
