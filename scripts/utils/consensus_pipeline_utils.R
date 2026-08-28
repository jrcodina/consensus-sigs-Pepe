library(dplyr)
library(rlang)

# Assign signature numbers
subgroup_to_sigs <- list(
  APOBEC                  = c("SBS2", "SBS13"),
  Clock_SBS1              = c("SBS1"),
  Clock_SBS5              = c("SBS5"),
  HRD                     = c("SBS3"),
  Tobacco                 = c("SBS4"),
  MMR                     = c("SBS6", "SBS14", "SBS15", "SBS20", "SBS21", "SBS26", "SBS44"),
  UV                      = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
  POLE                    = c("SBS10a", "SBS10b")
)

## Create Mesica Binary Tables 1 and 2

# Binary table 1 encompasses multiple signatures, while binary table 2 includes only the top signature

create_mesica_table <- function(signature_matrix, subgroup_mapping) {
  # Create an empty data frame with the SAMPLE_ID column
  mesica_table <- data.frame(SAMPLE_ID = signature_matrix$SAMPLE_ID)

  # Iterate through each subgroup and create the corresponding column in mesica_table
  for (subgroup in names(subgroup_mapping)) {
    # Extract the relevant signature names for the current subgroup
    sigs <- subgroup_mapping[[subgroup]]
    signature_matrix <- as.data.frame(signature_matrix)
    # Check if any of the signatures for the current subgroup have a value of 1 in each row
    mesica_table[[subgroup]] <- apply(signature_matrix[, sigs, drop = FALSE], 1, function(x) as.integer(any(x == 1)))
  }

  mesica_table
}
