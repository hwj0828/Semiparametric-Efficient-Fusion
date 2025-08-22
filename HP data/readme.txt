"HP_OBS_Control.csv": observational data that contains only control unit;
"HP_RCT_Control.csv": control unit in RCT data;
"HP_RCT_Treated.csv": treated unit in RCT data;
"hp_analysis.R": R code for data analysis.

Covariates: age, Pro.PRE, Education, Job, Gender, Gastritis.type, height, BMI, After.medical.treatement, Marriage, Nationality, Sick.or.vomit.PRE, Degree.of.stomach.pain.PRE, Hiccup.PRE, Acid.reflux.or.heartburn.PRE.

In the code:
Working model for propensity score is a logistic regression model by regression treatment on all covariates.
Working model for regression model is a logistic regression model by regressing outcome on all covariates. 
We use doubly robust estimator to estimate the average treatment effect when using only internal data.