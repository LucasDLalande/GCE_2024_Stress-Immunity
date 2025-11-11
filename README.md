# Early-life glucocorticoids accelerate lymphocyte count senescence in roe deer

Lucas D. Lalande<sup>1\*</sup>, Gilles Bourgoin<sup>1,2</sup>, Jeffrey Carbillet<sup>3</sup>, Louise Cheynel<sup>4</sup>, François Debias<sup>1</sup>, Hubert Ferté<sup>5</sup>, Jean-Michel Gaillard<sup>1</sup>, Rebecca Garcia<sup>1</sup>, Jean-François Lemaître<sup>1</sup>, Rupert Palme<sup>6</sup>, Maryline Pellerin<sup>7</sup>, Carole Peroz<sup>1,2</sup>, Benjamin Rey<sup>1</sup>, Pauline Vuarin<sup>1</sup>, Emmanuelle Gilot-Fromont<sup>1,2\*</sup>

<sup>1</sup> Université de Lyon, Université Lyon 1, CNRS, Laboratoire de Biométrie et Biologie Evolutive UMR 5558, F-69622 Villeurbanne, France\
<sup>2</sup> Université de Lyon, VetAgro Sup, 69280 Marcy l’Etoile, France\
<sup>3</sup> Current institution: Institute of Ecology and Earth Sciences, University of Tartu, 51014 Tartu, Estonia\
<sup>4</sup> Université de Lyon, Université Lyon 1, CNRS, Laboratoire d’Écologie des Hydrosystèmes Naturels et Anthropisés UMR 5023, F-69622 Villeurbanne, France\
<sup>5</sup> Université de Reims, Épidémio-Surveillance et Circulation de Parasites dans les Environnements UR 7510, 55100 Reims, France\
<sup>6</sup> Unit of Physiology, Pathophysiology and Experimental Endocrinology, Department of Biomedical Sciences, University of Veterinary Medicine, Veterinärplatz 1, 1210 Vienna, Austria\
<sup>7</sup> Office Français de la Biodiversité, Direction de la Recherche et de l’Appui Scientifique, Service Conservation et Gestion Durable des Espèces Exploités, 52210 Châteauvillain, France

Early-life glucocorticoids accelerate lymphocyte count senescence in roe deer\
*General and Comparative Endocrinology*, **2024**\
DOI: [10.1016/j.ygcen.2024.114595](https://doi.org/10.1016/j.ygcen.2024.114595)

<sup>\*</sup> Corresponding authors:\
[lalande.luke\@gmail.com](mailto:lalande.luke@gmail.com)\
[emmanuelle.gilotfromont\@vetagro-sup.fr](mailto:emmanuelle.gilotfromont@vetagro-sup.fr)

## Project description

We used these longitudinal data to test whether higher baseline GC levels measured within the first year of life were associated with a more pronounced immunosenescence and parasite susceptibility. We first assessed immunosenescence trajectories in these populations facing contrasting environmental conditions. Then, we found that juvenile GC levels can modulate lymphocyte trajectory. Lymphocyte depletion was accelerated late in life when FGMs were elevated early in life. Although the exact mechanism remains to be elucidated it could involve a role of GCs on thymic characteristics. In addition, elevated GC levels in juveniles were associated with a higher abundance of lung parasites during adulthood for individuals born during bad years, suggesting short-term negative effects of GCs on juvenile immunity, having in turn long-lasting consequences on adult parasite load, depending on juvenile environmental conditions. These findings offer promising research directions in assessing the carry-over consequences of GCs on life-history traits in the wild.

## Content

This deposit contains:

-   A description of the dataset used and where it can be found,
-   R scripts

to replicate the article mentioned above.

### Dataset

The current dataset has been obtained from a longitudinal survey as part of a Capture-Mark-Recapture program taking place in two french wild roe deer (*Capreolus capreolus*) populations since almost 50 years.

It can be found on Mendeley Data, DOI: [10.17632/bfjxsnj5h8.2](https://doi.org/10.17632/bfjxsnj5h8.2)

It includes:

-   idkit: the ID of the kit used for sampling
-   numind: Roe deer ID
-   pop: Population (TF: Trois-Fontaines, CH: Chizé)
-   cohorte: cohort
-   qualite_cohorte: cohort environmental quality defined as the mean fawn mass standardised for the date of capture
-   sexe: sex
-   annee: year of capture
-   datejulienne: Julian date of capture
-   qualite_an: year of capture environmental quality defined as the mean fawn mass standardised for the date of capture
-   ageannee: age (years)
-   issue_capt, lastyearcapt,lastobs,longevite: refer to the status of an individual (dead, alive, translocated to another population, …) and when it has been last observed or captured
-   masse: mass (kg)
-   id_bioch: ID for biochemical analyses
-   dureemin: time between capture and sampling (min)
-   WBC: white blood cells count
-   neutroP: Proportion of neutrophil
-   lymphoP: Proportion of lymphocyte
-   neutro: neutrophil count (10\^3 cells/mL)
-   eosino: eosinophil count (10\^3 cells/mL)
-   baso: basophil count (10\^3 cells/mL)
-   lympho: lymphocyte count (10\^3 cells/mL)
-   mono: monocyte count (10\^3 cells/mL)
-   alpha1: alpha1-globulin concentration (mg/mL)
-   alpha2: alpha2-globulin concentration (mg/mL)
-   beta: beta-globulin concentration (mg/mL)
-   gamma: gamma-globulin concentration (mg/mL)
-   HAP: haptoglobin concentration (mg/mL)
-   couleur_HAHL: coloration score of the serum
-   HA: haemagglutination score (10\^-2 log(dilution))
-   HL: haemolysis score (10\^-2 log(dilution))
-   SD: Gastrointestinal strongyles count
-   Trich: Trichuris species count
-   Coc: Coccidia count
-   Proto: Protostrongylid count
-   FGM_Priority: Whether faecal sampling is prioritize towards stress analyses or other analyses (nutrition, parasitism, …)
-   FGMngg: FGM levels (ng/g)

### R Scripts

-   `00_packages.R`: load all necessary libraries
-   `01_correction-FGM.R`: script for correction of FGM values according to the julian date of capture (linear and quadratic), to the time between capture and sampling, and to the time between sampling and freezing (immediate/\<24h)
-   `02_immunosenescence.R`: script to determine all senescence patterns (calls the `immunosenescence_function.R` file) and produce figures
-   `03_immunosenescence-FGM.R`: script to test for the effect of juvenile FGMs on retained senescence patterns and to produce figures
-   `immunosenescence_function.R`: function used in `02_immunosenescence.R` to determine patterns of senescence

## Configuration

Necessitate `R (>= 4.3.0)` and `RStudio`

## License

This repository contains the R code used to reproduce analyses from the article:

> Lalande et al. **(2024)**. Early-life glucocorticoids accelerate lymphocyte count senescence in roe deer. *General and Comparative Endocrinology.*

-   Code and scripts are released under the **CC-BY 4.0** license.
-   The data are available on Mendeley Data (DOI: [10.17632/bfjxsnj5h8.2](https://doi.org/10.17632/bfjxsnj5h8.2)) and are governed by their own license.

## Citation

To cite the article, please use :

> Lalande, Lucas D., Gilles Bourgoin, Jeffrey Carbillet, Louise Cheynel,François Debias, Hubert Ferté, Jean-Michel Gaillard, et al. **(2024)**.“Early-Life Glucocorticoids Accelerate Lymphocyte Count Senescencein Roe Deer.” *General and Comparative Endocrinology* 357: 114595. <https://doi.org/10.1016/j.ygcen.2024.114595>.

To cite the dataset and/or codes and scripts, please use :

> Lalande, Lucas D., Gilles Bourgoin, Jeffrey Carbillet, Louise Cheynel,François Debias, Hubert Ferté, Jean-Michel Gaillard, et al. **(2025)**.“Early-Life Glucocorticoids Accelerate Lymphocyte Count Senescencein Roe Deer.” *Mendeley Data*, V2. <https://doi.org/10.17632/BFJXSNJ5H8.1>.

## Contact:

For all queries :

Lucas D. Lalande - [lalande.luke\@gmail.com](mailto:lalande.luke@gmail.com)
