test_that("CompDbSource works", {
    expect_error(new("CompDbSource", dbfile = tempfile()), "unable")
    fl <- system.file("extdata", "MS1_example.txt", package = "MetaboAnnotation")
    expect_error(new("CompDbSource", dbfile = fl), "database")

    res <- new("CompDbSource")
    expect_true(validObject(res))

    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    res  <- new("CompDbSource", dbfile = fl)
    expect_true(validObject(res))

    expect_true(.validate_dbfile(fl))

    expect_output(show(res), "CompDbSource")
})

test_that("metadata,CompDbSource works", {
    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    src  <- new("CompDbSource", dbfile = fl)
    res <- metadata(src)
    expect_true(is.data.frame(res))
})

test_that("matchSpectra,Spectra,CompDbSource works", {
    fl <- system.file("sql", "CompDb.MassBank.sql", package = "CompoundDb")
    src  <- new("CompDbSource", dbfile = fl)

    res <- matchSpectra(pest_ms2, src, param = CompareSpectraParam())
    expect_s4_class(res, "MatchedSpectra")
    expect_equal(query(res), pest_ms2)
    expect_s4_class(target(res)@backend, "MsBackendDataFrame")
    expect_true(length(target(res)) == 0)

    library(CompoundDb)
    library(BiocParallel)
    register(SerialParam())
    if (file.exists(fl)) {
        qry <- Spectra(CompoundDb::CompDb(fl))[3]
        res <- matchSpectra(qry, src, param = CompareSpectraParam())
        expect_true(length(target(res)) == 4)
        expect_equal(MetaboAnnotation::matches(res)$target_idx, 1:4)
        expect_s4_class(target(res)@backend, "MsBackendDataFrame")
    }
})

test_that("MassBankSource works with AnnotationHub", {
    if (requireNamespace("AnnotationHub", quietly = TRUE)) {
        expect_error(MassBankSource(release = "other"), "not found")
        expect_error(MassBankSource(release = ""), "ambiguous")

        mb <- MassBankSource("2021.03")
        expect_s4_class(mb, "CompDbSource")
        expect_true(length(mb@dbfile) == 1L)
    }
})

test_that("MetFragSource works", {
    res <- MetFragSource()
    expect_s4_class(res, "MetFragSource")

    expect_error(MetFragSource("http://donotexist.com"), "404")
})

test_that("MetFragParam works", {
    res <- MetFragParam()
    expect_s4_class(res, "MetFragParam")
})

test_that(".peaks_to_metfrag_string works", {
    x <- cbind(mz = 1:4, intensity = 12:15)
    res <- .peaks_to_metfrag_string(x)
    expect_true(is.character(res))
    expect_equal(res, paste0(x[, 1L], "_", x[, 2L], collapse = ";"))

    x <- cbind(mz = numeric(), intensity = numeric())
    res <- .peaks_to_metfrag_string(x)
    expect_equal(res, character())
})

test_that("validAdducts works", {
    res <- validAdducts(MetFragParam())
    expect_true(is.list(res))
    expect_equal(names(res), c("positive", "negative"))

    DF <- DataFrame(rtime = c(1.3, 23.4))
    DF$mz <- list(c(12.1, 154, 155),
                  c(45.1, 456, 599))
    DF$intensity <- list(c(12, 54, 43),
                         c(1234, 454, 12334))
    sps <- Spectra(DF)
    expect_error(validAdducts(sps), "precursorAdduct")
    sps$precursorAdduct <- c("[M+H]+", "oooo")
    expect_error(validAdducts(sps), "polarity")
    sps$polarity <- c(1L, 0L)
    expect_error(validAdducts(sps), "adduct definitions")
    sps$precursorAdduct <- c("[M+H]+", "[M-H]-")
    expect_true(validAdducts(sps))
})

test_that(".param_to_list works", {
    res <- .param_to_list(MetFragParam())
    expect_true(is.list(res))
    expect_true(!any(names(res) %in% c("THRESHFUN", "validAdducts")))
})
