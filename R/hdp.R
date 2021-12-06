#' Hospital, doctor, patient (hdp) dataset
#'
#' This dataset has a three-level, hierarchical structure with patients nested within doctors within hospitals. The simulation code can be found at <https://stats.idre.ucla.edu/r/codefragments/mesimulation/#setup>.
#' @usage data(hdp)
#' @format A data frame with 8,525 rows and 27 variables:
#' \describe{
#'   \item{Age}{Continuous in years but recorded at a higher degree of accuracy.}
#'   \item{Married}{Binary, married/living with partner or single.}
#'   \item{FamilyHx}{Binary (yes/no), does the patient have a family history (Hx) of cancer?}
#'   \item{SmokingHx}{Categorical with three levels, current smoker, former smoker, never smoked.}
#'   \item{Sex}{Binary (female/male).}
#'   \item{CancerStage}{Categorical with four levels, stages 1-4.}
#'   \item{LengthofStay}{Count number of days patients stayed in the hospital after surgery.}
#'   \item{WBC}{Continuous, white blood count. Roughly 3,000 is low, 10,000 is middle, and 30,000 per microliter is high.}
#'   \item{RBC}{Continuous, red blood count..}
#'   \item{BMI}{Body mass index given by the formula (frac{kg}{meters^2}.}
#'   \item{IL6}{Continuous, interleukin 6, a proinflammatory cytokine commonly examined as an indicator of inflammation, cannot be lower than zero.}
#'   \item{CRP}{Continuous, C-reactive protein, a protein in the blood also used as an indicator of inflammation. It is also impacted by BMI.}
#'   \item{HID}{Hospital identifier.}
#'   \item{DID}{Doctor identifier}
#'   \item{Experience}{Years as a doctor.}
#'   \item{Lawsuits}{}
#'   \item{Medicaid}{}
#'   \item{School}{}
#'   \item{co2}{}
#'   \item{lungcapacity}{}
#'   \item{mobility}{}
#'   \item{nmorphine}{}
#'   \item{ntumors}{}
#'   \item{pain}{}
#'   \item{tumorsize}{}
#'   \item{wound}{}
#'   \item{remission}{Cancer in remission? 1 = yes, 0 = no.}
#' }
#' @source \url{https://stats.idre.ucla.edu/r/codefragments/mesimulation/#setup}
"hdp"
