x1 <- expand.grid(conds,conds)
x2 <- x1 %>%
  separate(Var1, sep = "_", into = c("Temp1", "Site1"), remove = F) %>%
  separate(Var2, sep = "_", into = c("Temp2", "Site2"), remove = F) %>%
  filter(Temp1 > Temp2 | Var1 == Var2)






##### create_conditions:  for comparison #####
create_conditions <- function(meta_data = deseq_data$meta_data){
  cat("Creating conditions for DEG comparison \n\n")
  conds <- unique(meta_data$condition)
  combinations <- combn(conds, 2,  simplify = F)
  conditions <- vector()
  count = 0
  for(comb in 1:length(combinations)){
    tmp1 <- strsplit(combinations[[comb]][1], split = "_")[[1]][2]
    tmp1.1 <- combinations[[comb]][1]
    tmp2 <- strsplit(combinations[[comb]][2], split = "_")[[1]][2]
    tmp2.1 <- combinations[[comb]][2]

    if(tmp1 == "AF" & tmp2 %in% c("ICN", "Protected","AF", "Exposed")){
      cat(tmp2, tmp1, "1\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp2 == "AF" & tmp1 %in% c("ICN", "Protected","AF", "Exposed")){
      cat(tmp1, tmp2, "2\n")
      name = paste("condition_", tmp1.1, "_vs_", tmp2.1, sep="")
    } else if(tmp1 == "Exposed" & tmp2 %in% c("Protected", "ICN")){
      cat(tmp2, tmp1, "3.1\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp2 == "Exposed" & tmp1 %in% c("Protected", "ICN", "Exposed")){
      cat(tmp1, tmp2, "3.2\n")
      name = paste("condition_", tmp1.1, "_vs_", tmp2.1, sep="")
    } else if(tmp2 == "Protected" & tmp1 %in% c("ICN", "Protected")){
      cat(tmp1, tmp2, "4\n")
      name = paste("condition_", tmp1.1, "_vs_", tmp2.1, sep="")
    } else if(tmp1 == "Protected" & tmp2 %in% c("ICN", "Protected")){
      cat(tmp2, tmp1, "4.2\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp1 == "ICN" & tmp2 %in% c("ICN")){
      cat(tmp2, tmp1, "5.2\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    } else if(tmp2 == "ICN" & tmp1 %in% c("ICN")){
      cat(tmp2, tmp1, "5.2\n")
      name = paste("condition_", tmp2.1, "_vs_", tmp1.1, sep="")
    }

    conditions <- append(name, conditions)
    count <- count + 1
  }
  cat("\t Returned ", length(conditions), " conditions \n")
  datatable(as.data.frame(conditions))
  return(conditions)
}
