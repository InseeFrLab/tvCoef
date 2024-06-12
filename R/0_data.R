#' Hansen Table
#' @docType data
#' @source Hansen, Bruce E. 1990. "Lagrange multiplier tests for parameter instability in non-linear models".
#' *University of Rochester*. <https://users.ssc.wisc.edu/~bhansen/papers/LMTests.pdf>.
"hansen_table"

#' French GDP
#'
#' Dataset containing the quarterly growth of the total gross domestic product (GDP)
#' of France and quarterly series of the French business climate in level and in difference.
#'
#' @details
#' Dataset containing the quarterly growth of the total gross domestic product (GDP, `"growth_gdp"`) of France,
#' in volumes chained at previous year prices, seasonally and working day adjusted;
#' and the French business climate in level.
#'
#' The French business climate is a monthly series, it is transformed into three quarterly series using the month's place in the quarter.
#' For example, `"bc_fr_m1"` contains the values in the first month of each quarter, and the `"diff_fr_m1"` is the difference
#' of the previous variable (the 2000Q1 value corresponds to the difference in business climate between January 2000 and October 1999).
#'
#' Data were downloaded March 15, 2024 and might therefore differ from the latest available data.
#' @docType data
#'
#' @format A quarterly \code{ts} object from 1949Q2 to 2024Q1.
#' @source INSEE
"gdp"

#' Business Surveys
#'
#' Dataset containing the quarterly growth of production in the manufacturing sector and its main sub-sectors,
#' quarterly balance of opinion of business surveys published by INSEE and Banque de France and
#' quarterly overhang of the industrial production index.
#'
#' @details
#'
#' Dataset containing the quarterly growth of production in the manufacturing sector and its main sub-sectors
#' and quarterly series series of business surveys published by INSEE and Banque de France.
#'
#' The sectors studied are:
#' - Manufacturing industry
#' - Food products and beverages (C1)
#' - Capital goods (C3)
#' - Transport equipments (C4)
#' - Other manufacturing (C5)
#'
#' `"manuf_prod"` contains the quarterly growth in production in the manufacturing sector,
#' and the sub-sectors are in the form `"prod_c1"`, `"prod_c3"`, `"prod_c4"` and `"prod_c5"`.
#'
#' The overhang of the industrial production index corresponds to the quarterly growth obtained extending the series by the last known value:
#' - `"overhang_ipi0"` is the quarterly growth obtained extending the series by the last value of the previous quarter (December, March, June, September);
#' - `"overhang_ipi1"` is the quarterly growth obtained extending the series by the first value of the current quarter (January, April, July, October);
#' - `"overhang_ipi2"` is the quarterly growth obtained extending the series by the second value of the current quarter (February, May, August, November).
#'
#' The business surveys being monthly, the balance of opinion are transformed into three quarterly series using the month's place in the quarter
#' (for example taking the values of January, April, July and October).
#' Variable names are constructed as the combination of several codes defined as follows:
#' - data source code (INSEE, `ins`, or Banque de France, `bdf`);
#' - name of the balance of opinion:
#'
#'   |**Code**          |**Definition**                  |
#'   |:-----------------|:-------------------------------|
#'   |bc                |Business climate                |
#'   |oscd              |Overall order books             |
#'   |tppa and prodpas  |Past production                 |
#'   |tppre and prodpre |Personal production expectation |
#'   |sitcar            |Situation of order books        |
#'   |evocar            |Evolution of order books        |
#'   |prix              |Selling prices                  |
#'   |stocks            |Inventories of finished goods   |
#'   |tres              |Cash position                   |
#'   |tuc               |Capacity utilisation rate       |
#'
#' - sector: nothing for the manufacturing industri and `"c1"`, `"c3"`, `"c4"` or `"c5"` for the sub-sectors;
#' - place of the month in the quarter: `m1`, `m2` or `m3` for the first, second or third month of the quarter.
#'
#' The dataset also contains some dummies labelled `"indYYYYQX"`, where `YYYY` is the year and `X` is the quarter.
#'
#' @docType data
#' @format A quarterly \code{ts} object from 1949Q2 to 2024Q1.
#'
#' @source INSEE, Banque de France
"manufacturing"
