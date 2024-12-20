###################################################
##### model_benchmark.r
##### Last Modified by S.Zhang on 2024/11/24
###################################################

library(dplyr)
library(data.table)
library(caret)
library(pROC)
library(randomForest)

p_data <- "../data"
p_data_mimic_iv <- "../data/mimic_iv"


#####----------------------------------------------------------------------#####
##### Data Processing #####
#####----------------------------------------------------------------------#####

##### Load and prepare Data #####
df_analytical_ready <- read.csv(file.path(p_data, "analytical_ready.csv"), stringsAsFactors = F)
names(df_analytical_ready)

# keep only records for sepsis patients only before onset
df_analytical_ready <- df_analytical_ready%>%ungroup()%>%
  mutate(sepsis_char = as.character(sepsis))



# Split the data into training (80%) and test (20%) sets
set.seed(123) # For reproducibility
df_patient_unique <- df_analytical_ready%>%
  select(unique_id, sepsis)%>%unique()
trainIndex <- createDataPartition(df_patient_unique$sepsis, p = 0.80, list = FALSE)
patientid_train <- df_patient_unique[trainIndex,]
patientid_test <- df_patient_unique[-trainIndex,]
trainData <- df_analytical_ready%>%inner_join(patientid_train)
testData <- df_analytical_ready%>%inner_join(patientid_test)
table(trainData$sepsis)
table(testData$sepsis)


#####----------------------------------------------------------------------#####
##### qSOFA #####
#####----------------------------------------------------------------------#####
testData_qsofa <- testData%>%
  mutate(flag_qsofa = ifelse(resprate_r >= 22 & sbp_r <= 100, 1, 0),
         flag_qsofa = ifelse(is.na(flag_qsofa), 0, flag_qsofa))%>%
  group_by(unique_id)%>%
  summarise(sepsis = max(sepsis),
            flag_qsofa = max(flag_qsofa))

table(testData_qsofa$sepsis, testData_qsofa$flag_qsofa)

accuracy_qsofa <- (800 + 23) / (21 + 55 + 800 +23)
recall_qsofa <- 23/(23+21)

#####----------------------------------------------------------------------#####
##### Apply logistic regression with interaction term #####
#####----------------------------------------------------------------------#####

model <- glm(as.factor(sepsis) ~ gender + race_r + 
              time_encoding_1 + time_encoding_2 + time_encoding_3 + time_encoding_4 + 
              time_encoding_5 + time_encoding_6 + time_encoding_7 + time_encoding_8 + 
              temperature_1 + temperature_2 + temperature_3 + temperature_4 + 
              temperature_5 + temperature_6 + temperature_7 + temperature_8 + 
              heartrate_1 + heartrate_2 + heartrate_3 + heartrate_4 + 
              heartrate_5 + heartrate_6 + heartrate_7 + heartrate_8 + 
              resprate_1 + resprate_2 + resprate_3 + resprate_4 + resprate_5 + resprate_6 + resprate_7 + resprate_8 + 
              o2sat_1 + o2sat_2 + o2sat_3 + o2sat_4 + o2sat_5 + o2sat_6 + o2sat_7 + o2sat_8 +    
              sbp_1 + sbp_2 + sbp_3 + sbp_4 + sbp_5 + sbp_6 + sbp_7 + sbp_8 + 
              dbp_1 + dbp_2 + dbp_3 + dbp_4 + dbp_5 + dbp_6 + dbp_7 + dbp_8, 
               data = trainData, family = binomial)


predicted_probabilities <- predict(model, newdata = testData, type = "response")
predicted_outcome <- ifelse(predicted_probabilities > 0.5, 1, 0)

testData$predicted <- predicted_outcome

# Summarize the predicted data to patient level
testData_patient <- testData%>%
  filter(!is.na(predicted))%>%
  group_by(unique_id)%>%
  mutate(predicted_patient = max(predicted, na.rm = T))%>%
  select(unique_id, sepsis, predicted_patient)%>%
  unique()


# Calculate the AUC
roc_obj <- roc(testData_patient$sepsis, testData_patient$predicted_patient)
auc_value <- auc(roc_obj)
auc_value

# evaluation
confusion_matrix <- table(Predicted = testData_patient$predicted_patient, Actual = testData_patient$sepsis)
confusion_matrix


# Calculate accuracy
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(paste("Accuracy:", accuracy))
# Calculate precision
precision <- confusion_matrix[2, 2] / sum(confusion_matrix[2, ])
print(paste("Precision:", precision))
# Calculate recall (sensitivity)
recall <- confusion_matrix[2, 2] / sum(confusion_matrix[, 2])
print(paste("Recall:", recall))

# Calculate F1-score
f1_score <- 2 * (precision * recall) / (precision + recall)
print(paste("F1-Score:", f1_score))



#####----------------------------------------------------------------------#####
##### Random Forest
#####----------------------------------------------------------------------#####

trainData_clean <- trainData%>%
  mutate(sepsis_factor = as.factor(sepsis_char))%>%
  select(sepsis_factor, gender, race_r, 
           time_encoding_1, time_encoding_2, time_encoding_3, time_encoding_4, 
           time_encoding_5, time_encoding_6, time_encoding_7, time_encoding_8, 
           temperature_1, temperature_2, temperature_3, temperature_4, 
           temperature_5, temperature_6, temperature_7, temperature_8, 
           heartrate_1, heartrate_2, heartrate_3, heartrate_4, 
           heartrate_5, heartrate_6, heartrate_7, heartrate_8, 
           resprate_1, resprate_2, resprate_3, resprate_4, resprate_5, resprate_6, resprate_7, resprate_8, 
           o2sat_1, o2sat_2, o2sat_3, o2sat_4, o2sat_5, o2sat_6, o2sat_7, o2sat_8,    
           sbp_1, sbp_2, sbp_3, sbp_4, sbp_5, sbp_6, sbp_7, sbp_8, 
           dbp_1, dbp_2, dbp_3, dbp_4, dbp_5, dbp_6, dbp_7, dbp_8)
trainData_clean <- na.omit(trainData_clean)

testData_clean <- testData%>%
  mutate(sepsis_factor = as.factor(sepsis_char))%>%
  select(sepsis_factor, gender, race_r, 
         time_encoding_1, time_encoding_2, time_encoding_3, time_encoding_4, 
         time_encoding_5, time_encoding_6, time_encoding_7, time_encoding_8, 
         temperature_1, temperature_2, temperature_3, temperature_4, 
         temperature_5, temperature_6, temperature_7, temperature_8, 
         heartrate_1, heartrate_2, heartrate_3, heartrate_4, 
         heartrate_5, heartrate_6, heartrate_7, heartrate_8, 
         resprate_1, resprate_2, resprate_3, resprate_4, resprate_5, resprate_6, resprate_7, resprate_8, 
         o2sat_1, o2sat_2, o2sat_3, o2sat_4, o2sat_5, o2sat_6, o2sat_7, o2sat_8,    
         sbp_1, sbp_2, sbp_3, sbp_4, sbp_5, sbp_6, sbp_7, sbp_8, 
         dbp_1, dbp_2, dbp_3, dbp_4, dbp_5, dbp_6, dbp_7, dbp_8)
testData_clean <- na.omit(testData_clean)


# Train a random forest model
set.seed(123)  # For reproducibility
rf_model <- randomForest(sepsis_factor ~ ., data=trainData_clean, ntree=100, mtry=2)

# Print model summary
print(rf_model)

# Predict on test set
predicted_probs <- predict(rf_model, newdata = testData_clean, type = "prob")

df_predicted_probs <- data.frame(predicted_probs)
df_predicted_probs$actual <- testData_clean$sepsis_factor
df_predicted_probs <- df_predicted_probs%>%
  mutate(predict_5 = ifelse(X1>0.5, 1, 0))
table(Predicted = df_predicted_probs$predict_5, Actual = df_predicted_probs$actual)


# Add predicted value back to the original dataset
testData$predict_sepsis <- df_predicted_probs$predict_5

testData_clean_patients <- testData%>%
  group_by(unique_id)%>%
  summarise(sepsis_actual = max(as.numeric(sepsis)),
            sepsis_predict = max(predict_sepsis))



# Assuming "1" is the positive class
roc_curve <- roc(testData_clean_patients$sepsis_actual, testData_clean_patients$sepsis_predict)

# Plot the ROC curve
plot(roc_curve, main = "ROC Curve for Random Forest Model")
auc_value <- auc(roc_curve)
print(paste("AUC:", auc_value))


# evaluation
confusion_matrix_rf <- table(Predicted = testData_clean_patients$sepsis_predict, Actual = testData_clean_patients$sepsis_actual)
confusion_matrix_rf


# Calculate accuracy
accuracy <- sum(diag(confusion_matrix_rf)) / sum(confusion_matrix_rf)
print(paste("Accuracy:", accuracy))
# Calculate precision
precision <- confusion_matrix_rf[2, 2] / sum(confusion_matrix_rf[2, ])
print(paste("Precision:", precision))
# Calculate recall (sensitivity)
recall <- confusion_matrix_rf[2, 2] / sum(confusion_matrix_rf[, 2])
print(paste("Recall:", recall))


