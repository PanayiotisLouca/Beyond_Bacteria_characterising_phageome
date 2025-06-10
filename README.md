# Beyond Bacteria: Characterising the Gut Phageome in the General Population 
## Authors 
*Panayiotis Louca, Mohammadali Khan Mirzaei, Afroditi Kouraki, Xue Peng, Erfan Khamespanah, Yu Lin, Robert Pope, Alessia Visconti, Francesco Asnicar, Daniel Kirk, Ricardo Costeira, Nicola Segata, Mario Falchi, Jordana T Bell, Tim D. Spector, Lindsey A Edwards, Li Deng, Ana M. Valdes, Cristina Menni*
  
---
  This repository contains the R code and analysis scripts used in our comprehensive study of the human gut phageome in the general population (DOI: To be updated). 
---

  ## ℹ️ Repository Description 
  
  - **Phage-Taxa Network Analysis**
      - `Phage_taxa_ggraph_network_strong_pos_corr_SCRIPT.R`: Constructs network for strong positive correlations (Figure 3A).
      - `Phage_taxa_ggraph_network_mod_neg_corr_SCRIPT.R`: Builds network for moderate negative correlations between phage taxa (Figure 3B).
  - **Phage-Microbial Metabolite Associations**
      - `Phage_SCFA_CREATE_SCRIPT.R`: Tests associations between viral contigs & SCFAs in serum and stool using linear mixed effect models adjusting for age, sex, and BMI as fixed effect variables, and family relatedness and batch as random effects.
      - `Phage_bile_acids_CREATE_SCRIPT.R`: Tests associations between viral contigs & secondary bile acids from serum and stool using linear mixed effect models adjusting for age, sex, and BMI as fixed effect variables, and family relatedness and batch as random effects.
  - **Phage-Metabolite Network Construction**
      - `Phage_metabolite_ggraph_network_SCRIPT.R`: Constructs a correlation network between signficiantly associated viral contigs and serum & stool SCFAs & secondary bile acids.
  - **Phage-Diet Associations**
      - `Phage_diet_CREATE_SCRIPT.R`: Tests associations between viral contigs & alpha diveristy of the gut phageome and dietary nutrient intakes and the healthy eating index using linear mixed effect models. 
  - **Phage-Metabolic Health Associations**
      - `Phage_metabolic_health_CREATE_SCRIPT.R`: Tests associations between viral contigs & metabolic health parameters (Triglycerides, TyG index, & glucose) using linear mixed effect modelling.
  

  ## 📂 Repository Structure
```
├── Phage_Taxa_Network_Analysis
│   ├── Phage_taxa_ggraph_network_strong_pos_corr_SCRIPT.R
│   ├── Phage_taxa_ggraph_network_mod_neg_corr_SCRIPT.R
├── Phage_Microbial_Metabolite_Associations
│   ├── Phage_SCFA_CREATE_SCRIPT.R
│   ├── Phage_bile_acids_CREATE_SCRIPT.R
├── Phage_Metabolite_Network_Construction
│   ├── Phage_metabolite_ggraph_network_SCRIPT.R
├── Phage_Diet_Associations
│   ├── Phage_diet_CREATE_SCRIPT.R
└── Phage_Metabolic_Health_Associations
    └── Phage_metabolic_health_CREATE_SCRIPT.R
```
          
  ## Citation 
  If you use this code, please cite:
  Louca, P. et al. (2025). Beyond Bacteria: Characterising the Gut Phageome in the General Population. [To be updated].
  
  
  
