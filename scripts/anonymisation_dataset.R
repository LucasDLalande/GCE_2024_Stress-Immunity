data <- read.csv("dataset/00_dataset-immunosenescence-FGM.csv", h=T, sep=";", dec=".")

# Numind

data_unique_ID <- unique(factor(data$numind))

anonymID_correspondance <- data.frame("numind"=unique(factor(data$numind)), "anonym_numind"=paste0("ID", seq_along(data_unique_ID)))

write.table(anonymID_correspondance, "dataset/dataset_anonym/anonymID_correspondance.csv", row.names = F, sep=";", dec=".")

data_anonymID <- merge(data, anonymID_correspondance, by="numind", all=T)

# IDkit

data_unique_IDkit <- unique(factor(data$idkit))

anonymIDkit_correspondance <- data.frame("idkit"=unique(factor(data$idkit)), "anonym_idkit"=paste0("IDkit", seq_along(data_unique_IDkit)))

write.table(anonymIDkit_correspondance, "dataset/dataset_anonym/anonymIDkit_correspondance.csv", row.names = F, sep=";", dec=".")

data_anonymID_IDkit <- merge(data_anonymID, anonymIDkit_correspondance, by="idkit", all=T)

# ID_bioch

data_unique_IDbioch <- unique(factor(data$id_bioch))

anonymIDbioch_correspondance <- data.frame("id_bioch"=unique(factor(data$id_bioch)), "anonym_id_bioch"=paste0("ID_bioch", seq_along(data_unique_IDbioch)))

write.table(anonymIDbioch_correspondance, "dataset/dataset_anonym/anonymIDbioch_correspondance.csv", row.names = F, sep=";", dec=".")

data_anonymID_IDkit_IDbioch <- merge(data_anonymID_IDkit, anonymIDbioch_correspondance, by="id_bioch", all=T)



write.table(data_anonymID_IDkit_IDbioch, "dataset/dataset_anonym/00_dataset-immunosenescence-FGM_anonym.csv", row.names = F, sep=";", dec=".")
