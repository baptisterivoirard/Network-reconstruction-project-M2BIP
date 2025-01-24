#librairies

library(miic)
library(igraph)

#Dataset METABRIC utilisé : https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric

# ouverture du fichier tableur

data=read.csv(file.choose(),header=TRUE,sep=',',row.names = NULL)

#création du dataframe final pour le premier réseau

data_test=data.frame(data$age_at_diagnosis,data$cancer_type_detailed,
                     data$er_status,data$neoplasm_histologic_grade,
                     data$inferred_menopausal_state,
                     
                     data$overall_survival_months,data$overall_survival,
                     data$pr_status,data$tumor_size,
                     data$tumor_stage,data$death_from_cancer,
                     
                     data$cdkn1b,data$cdkn2a,
                     data$foxo1,data$foxo3,data$foxp1,
                     data$tp53,
                     data$brca1,data$brca2,
                     data$akt1,data$akt2,
                     data$pik3ca,
                     data$ahnak2,data$ahnak,
                     data$syne1,
                     data$kmt2d,
                     data$gata3,
                     
                     data$rb1,
                     data$myc,
                     
                     data$cdk1,data$cdk2,
                     data$cdk4,data$cdk6,data$cdk8)

# écriture du data frame dans un fichier tableur

write.csv(data_test,"C:/Users/Paul/Desktop/M2/MU5IN754_RESYS/data_reseau_1_f.csv", row.names=FALSE)



## Calcul de la mutual information et création de fichier tableur pour la reconstruction des sous réseaux :

# Exemple pour tumor size : Sélection des gènes avec la plus grande mutual information avec tumor size : 
# Séléction de toute les data d'expression de gène
subdata_tumorsize = + data[, 32: 520]
subdata_tumorsize$tumor_size <- data$tumor_size

# Initialisation d'une liste pour stocker les résultats d'information mutuelle
mutual_info_results <- list()

# Parcourir chaque colonne de gène dans 'subdata_survival' (en excluant 'overall_survival')
for (gene in names(subdata_tumorsize)[names(subdata_tumorsize) != "tumor_size"]) {
  
  # Calculer l'information mutuelle entre chaque gène et la variable clinique 'overall_survival'
  res <- computeMutualInfo(
    x = subdata_tumorsize[[gene]], 
    y = subdata_tumorsize$tumor_size, 
    maxbins = 10,                cplx = "nml",             
    plot = FALSE              )
  
  # Stocker le nom du gène et l'information mutuelle
  mutual_info_results[[gene]] <- res$info
}

# Transformer les résultats en dataframe pour les trier
mi_df <- data.frame(
  gene = names(mutual_info_results),
  mutual_info = unlist(mutual_info_results)
)

# Classer les gènes par ordre décroissant d'information mutuelle
mi_df <- mi_df[order(-mi_df$mutual_info), ]

# Afficher les résultats
print(mi_df)

# Séléction uniquement des gènes avec MI > 0
genes_with_positive_mi <- mi_df$gene[mi_df$mutual_info > 0]

data_test <- data.frame(data[, genes_with_positive_mi],data$tumor_size)

# Crée le fichier .csv pour la reconstruction du réseau par le serveur MIIC
write.csv(data_test, 'C:\\Users\\rivoi\\Documents\\M2MAD\\REYS\\Project\\Project.csv', row.names = FALSE)


