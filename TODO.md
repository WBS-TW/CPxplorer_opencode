\# TODO



Already done:
* For the Advanced settings, refactor the code for user input of transformation products to make it easier and more systematic to add more transformation products. Suggest additional transformation products that are relevant besides the current ones. More examples are in the @/inst/CPions_TP_formula.xlsx file. When you add additional transformation products, also fill in these in the rows of the @inst/CPions_TP_formula.xlsx file.
Add a checkbox (unchecked by default) which when the user check it and will, instead of prefixed selection, add a textInput where the user put in one or several formula for transformation products. The grammar for the transformation product insertion is recommended to be in the same style as for the predefined transformation product: -H+OH; -Cl+OH; -2H+2OH; where -H refers to one H atom being replaced and +OH is the replacement with an OH. The ";" sign is a separator to mark addition of a new transformation product. A numeric prefix before the atom such as -2H indicate loss of 2 H atoms, and +2OH indicate addition of 2 OH molecules. You must consideration the formula feasibility of the user input, for example, -2H+O, is feasible which indicate replacement of two hydrogens to the same carbon atom with a double bond to oxygen.

The user input should be valid for the basic formula of PCA, PCO and BCA. Verify the results from the chemical formula with the different transformation products that these are valid.

