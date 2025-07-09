import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import LabelEncoder
# Removed ML imports as requested: RandomForestClassifier, classification_report, confusion_matrix
from sklearn.cluster import KMeans, DBSCAN # Added DBSCAN
import plotly.express as px
import warnings
from datetime import datetime
import os
import base64 # Import base64 for the download button
import tempfile
import shutil # Import shutil for directory operations
import statsmodels.api as sm

import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
from sklearn.metrics import roc_curve, roc_auc_score

# Set up the app
st.set_page_config(layout="wide", page_title="Respiratory Virus Analysis")
warnings.filterwarnings('ignore')
pd.set_option('display.max_columns', None)
plt.style.use('ggplot')

# Custom CSS for better styling
st.markdown("""
    <style>
        .main {padding: 2rem;}
        .stButton>button {width: 100%;}
        .stDownloadButton>button {width: 100%;}
        .stFileUploader>div>div>div>button {width: 100%;}
        .reportview-container .main .block-container {padding-top: 2rem; padding-bottom: 2rem;}
        h1 {color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;}
        h2 {color: #2980b9; margin-top: 30px;}
        h3 {color: #16a085;}
        .sidebar .sidebar-content {background-color: #f8f9fa;}
    </style>
""", unsafe_allow_html=True)

# Define the path to your Excel file
# Make sure 'qPCR_wide_by_agent.xlsx' is in the same directory as your Streamlit script
# Or provide a full path like r"C:\path\to\your_dataset.xlsx"
DATA_FILE_PATH = 'Ana_Veri_Guncellenmis_Son_Hali.xlsx'

# Sidebar for analysis options
with st.sidebar:
    st.title("Analysis Settings")

    st.success(f"Using pre-loaded file: {DATA_FILE_PATH}") # Inform the user

    # Updated analysis options: removed ML and Prediction, kept Advanced Statistical Analysis
    analysis_options = st.multiselect(
        "Select analysis sections to include:",
        ["Demographic Analysis", "Virus Distribution", "Time Series Analysis",
         "Environmental Factors", "Co-infections", "Advanced Statistical Analysis", "Clustering"],
        default=["Demographic Analysis", "Virus Distribution", "Time Series Analysis", "Advanced Statistical Analysis"]
    )

    show_raw_data = st.checkbox("Show raw data preview")
    show_advanced = st.checkbox("Show advanced options")

    if show_advanced:
        st.subheader("Advanced Options")
        min_support = st.slider("Minimum support for association rules", 0.001, 0.1, 0.005, 0.001)
        # random_state for clustering is now handled within the clustering section
        # n_clusters for K-Means is now handled within the clustering section

# Main app content
st.title("Respiratory Virus Analysis Dashboard")

# Load and process data
@st.cache_data
def load_data(file_path):
    """
    Loads and preprocesses the respiratory virus data from an Excel file.

    Args:
        file_path (str): The path to the Excel data file.

    Returns:
        tuple: A tuple containing:
            - df (pd.DataFrame): The processed DataFrame.
            - actual_virus_columns (list): List of identified virus column names.
            - env_temp_col (str): Name of the temperature column.
            - env_hum_col (str): Name of the humidity column.
            - labels (list): List of age group labels.
    """
    try:
        df = pd.read_excel(file_path)

        # Standardize column names by stripping whitespace, replacing spaces with underscores,
        # and replacing slashes with underscores. Also remove parentheses.
        df.columns = df.columns.str.strip().str.replace(' ', '_').str.replace('/', '_').str.replace('(', '').str.replace(')', '')


        # --- MODIFICATION: Update date column identification to include 'Çalisma_Tar.' ---
        date_cols_priority = ['Çalisma_Tar.', 'Uygulama_Tarihi', 'Study_Date', 'Sample_Date', 'Tarih']
        actual_date_col = None
        for col in date_cols_priority:
            if col in df.columns:
                df[col] = pd.to_datetime(df[col], errors='coerce')
                # Check if there are any non-NA values after conversion
                if df[col].notna().any():
                    actual_date_col = col
                    break

        if actual_date_col:
            # Drop rows where the identified date column is NaN
            df.dropna(subset=[actual_date_col], inplace=True)
            df['Analyzed_Date'] = df[actual_date_col]
            df['Year'] = df['Analyzed_Date'].dt.year
            df['Month'] = df['Analyzed_Date'].dt.month
            # Prioritize 'Hafta' column if it exists and is numeric, otherwise calculate from Analyzed_Date
            if 'Hafta' in df.columns:
                # Convert 'Hafta' to numeric, coerce errors to NaN, then fill NaNs
                # with week calculated from 'Analyzed_Date'
                df['Week_of_Year'] = pd.to_numeric(df['Hafta'], errors='coerce').fillna(df['Analyzed_Date'].dt.isocalendar().week.astype(int))
            else:
                df['Week_of_Year'] = df['Analyzed_Date'].dt.isocalendar().week.astype(int)
        else: # Fallback if no specific date column is found
            if 'Year' in df.columns and 'Month' in df.columns: # Using 'Year' and 'Month' as per new headers
                try:
                    # Attempt to create 'Analyzed_Date' from 'Year' and 'Month'
                    df['Analyzed_Date'] = pd.to_datetime(df['Year'].astype(str) + '-' + df['Month'].astype(str) + '-01', errors='coerce')
                    df.dropna(subset=['Analyzed_Date'], inplace=True)
                    df['Year'] = df['Analyzed_Date'].dt.year
                    df['Month'] = df['Analyzed_Date'].dt.month
                    if 'Hafta' in df.columns: # Use 'Hafta' if available
                        df['Week_of_Year'] = pd.to_numeric(df['Hafta'], errors='coerce').fillna(df['Analyzed_Date'].dt.isocalendar().week.astype(int))
                    else:
                        df['Week_of_Year'] = df['Analyzed_Date'].dt.isocalendar().week.astype(int)
                except Exception:
                    st.warning("Could not create 'Analyzed_Date' from 'Year' and 'Month'. Time series analysis might be limited.")
                    # Assign NaT/NaN if date creation fails to prevent further errors
                    df['Analyzed_Date'] = pd.NaT
                    df['Year'] = np.nan
                    df['Month'] = np.nan
                    df['Week_of_Year'] = np.nan
            else:
                st.warning("No suitable date column ('Çalisma_Tar.', 'Uygulama_Tarihi', 'Study_Date', 'Sample_Date', 'Tarih', or 'Year'/'Month') found. Time series analysis will be skipped.")
                # Assign NaT/NaN if no date columns are found
                df['Analyzed_Date'] = pd.NaT
                df['Year'] = np.nan
                df['Month'] = np.nan
                df['Week_of_Year'] = np.nan

        # --- NEW: Create Season column ---
        if 'Month' in df.columns:
            def get_season(month):
                if month in [12, 1, 2]:
                    return 'Winter'
                elif month in [3, 4, 5]:
                    return 'Spring'
                elif month in [6, 7, 8]:
                    return 'Summer'
                elif month in [9, 10, 11]:
                    return 'Autumn'
                return 'Unknown'
            df['Season'] = df['Month'].apply(get_season)
        else:
            st.warning("Month column not found, Season column will not be created.")
            df['Season'] = 'Unknown'


        # Process age
        labels = ['0-5', '6-18', '19-40', '41-65', '65+'] # Define labels for age groups
        if 'Yaş' in df.columns:
            df['Yaş'] = pd.to_numeric(df['Yaş'], errors='coerce') # Convert 'Yaş' to numeric, coercing errors
            df.dropna(subset=['Yaş'], inplace=True) # Drop rows with NaN in 'Yaş'
            bins = [0, 5, 18, 40, 65, 120] # Define age bins
            df['Yaş_Grubu'] = pd.cut(df['Yaş'], bins=bins, labels=labels, right=False) # Create age groups
        else:
            st.warning("Age column ('Yaş') not found. Age-based analysis will be skipped.")
            df['Yaş_Grubu'] = 'Unknown' # Placeholder if 'Yaş' column is missing

        # --- MODIFICATION: Update environmental factor column names to match provided headers ---
        env_temp_col = 'Temp_ay_ortalaması' # Use 'Temp_ay_ortalaması'
        env_hum_col = 'Humidity_ayortalaması' # Use 'Humidity_ayortalaması'

        if env_temp_col in df.columns:
            df[env_temp_col] = pd.to_numeric(df[env_temp_col], errors='coerce')
            df.dropna(subset=[env_temp_col], inplace=True)
        else:
            df[env_temp_col] = np.nan # Ensure column exists even if not found
            st.warning(f"Temperature column ('{env_temp_col}') not found. Environmental analysis might be limited.")

        if env_hum_col in df.columns:
            df[env_hum_col] = pd.to_numeric(df[env_hum_col], errors='coerce')
            df.dropna(subset=[env_hum_col], inplace=True)
        else:
            df[env_hum_col] = np.nan # Ensure column exists even if not found
            st.warning(f"Humidity column ('{env_hum_col}') not found. Environmental analysis might be limited.")

        # Identify virus columns from a predefined template list
        # --- MODIFICATION: Update virus column names to match provided headers ---
        virus_columns_template = [
            'Adenovirus', 'Coronavirus_HKU1', 'Enterovirus_Solunum',
            'Human_Bocavirus', 'Human_Coronavirus_229E', 'Human_Coronavirus_NL63',
            'Human_Coronavirus_OC43', 'Human_Metapneumovirus',
            'Human_parechovirus_Solunum', 'Influenza_A', 'Influenza_B',
            'Parainfluenza_Virus_1', 'Parainfluenza_Virus_2',
            'Parainfluenza_Virus_3', 'Parainfluenza_Virus_4',
            'Respiratuvar_sinsityal_virüs_A_B', 'Rhinovirus', 'SARS-COV-2'
        ]

        actual_virus_columns = [col for col in virus_columns_template if col in df.columns]

        if actual_virus_columns:
            # Calculate 'Total_Virus_Count' by summing positive virus detections for each row
            df['Total_Virus_Count'] = df[actual_virus_columns].sum(axis=1)
        else:
            st.error("No virus columns found in the dataset. Virus-based analysis will be skipped.")
            df['Total_Virus_Count'] = 0 # Placeholder if no virus columns are found

        # Final cleanup for analysis sections: drop rows with NaNs in critical columns
        cols_to_check = []
        if 'Yaş_Grubu' in df.columns:
            cols_to_check.append('Yaş_Grubu')
        cols_to_check.extend([col for col in actual_virus_columns if col in df.columns])

        if cols_to_check: # Only apply if there are columns to check
             df.dropna(subset=cols_to_check, inplace=True)

        if len(df) < 10:
            st.warning(f"After preprocessing, only {len(df)} records remain. Some analyses might be affected.")

        return df, actual_virus_columns, env_temp_col, env_hum_col, labels

    except FileNotFoundError:
        st.error(f"Error: The data file '{file_path}' was not found. Please ensure it's in the correct directory.")
        st.stop()
    except Exception as e:
        st.error(f"Error loading or processing data: {e}")
        st.stop()

# Call load_data with the predefined file path
df, actual_virus_columns, env_temp_col, env_hum_col, labels = load_data(DATA_FILE_PATH)

# --- NEW CHECK: Stop if DataFrame is empty after loading ---
if df.empty:
    st.error("The DataFrame is empty after data loading and preprocessing. Please check your Excel file and the column names. If the issue persists, consider checking for hidden characters in column names or data types.")
    st.stop()
# --- END NEW CHECK ---

if show_raw_data:
    st.subheader("Raw Data Preview")
    st.dataframe(df.head())

# Analysis sections (existing sections remain largely unchanged, just ensuring column names are correct)
if "Demographic Analysis" in analysis_options:
    st.header("1. Demographic Analysis")

    col1, col2 = st.columns(2)

    with col1:
        if 'Cinsiyet' in df.columns:
            st.subheader("Gender Distribution")
            fig, ax = plt.subplots(figsize=(8, 6))
            df['Cinsiyet'].value_counts().plot(kind='pie', autopct='%1.1f%%',
                                              colors=['skyblue', 'lightcoral'], ax=ax)
            ax.set_ylabel('')
            st.pyplot(fig)
        else:
            st.warning("Gender column ('Cinsiyet') not found.")

    with col2:
        if 'Yaş' in df.columns:
            st.subheader("Age Distribution")
            fig, ax = plt.subplots(figsize=(8, 6))
            sns.histplot(df['Yaş'], bins=30, kde=True, color='teal', ax=ax)
            ax.set_xlabel('Age')
            ax.set_ylabel('Patient Count')
            st.pyplot(fig)
        else:
            st.warning("Age column ('Yaş') not found.")

    if 'Yaş_Grubu' in df.columns and df['Yaş_Grubu'].nunique() > 1: # Ensure there's more than one age group
        st.subheader("Age Group Distribution")
        fig, ax = plt.subplots(figsize=(10, 6))
        df['Yaş_Grubu'].value_counts().sort_index().plot(kind='bar', color='mediumseagreen', ax=ax)
        ax.set_xlabel('Age Group')
        ax.set_ylabel('Patient Count')
        ax.tick_params(axis='x', rotation=45)
        st.pyplot(fig)
    else:
        st.warning("Age Group column ('Yaş_Grubu') not found or only one group present.")

if "Virus Distribution" in analysis_options and actual_virus_columns:
    st.header("2. Virus Distribution Analysis")

    st.subheader("Virus Case Counts")
    virus_counts = df[actual_virus_columns].sum().sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(12, 8))
    virus_counts.plot(kind='barh', color='darkorange', ax=ax)
    ax.set_xlabel('Case Count')
    ax.set_ylabel('Virus Type')
    st.pyplot(fig)

    col1, col2 = st.columns(2)

    with col1:
        if 'Yaş_Grubu' in df.columns and df['Yaş_Grubu'].nunique() > 1:
            st.subheader("Virus Distribution by Age Group")
            fig, ax = plt.subplots(figsize=(10, 8))
            age_virus = df.groupby('Yaş_Grubu')[actual_virus_columns].sum().T
            sns.heatmap(age_virus, annot=True, fmt='d', cmap='YlOrRd', ax=ax)
            ax.set_xlabel('Age Group')
            ax.set_ylabel('Virus Type')
            st.pyplot(fig)
        else:
            st.warning("Age Group column ('Yaş_Grubu') not found or only one group present for virus distribution.")

    with col2:
        if 'Cinsiyet' in df.columns and df['Cinsiyet'].nunique() > 1:
            st.subheader("Virus Distribution by Gender")
            fig, ax = plt.subplots(figsize=(10, 8))
            gender_virus = df.groupby('Cinsiyet')[actual_virus_columns].sum().T
            sns.heatmap(gender_virus, annot=True, fmt='d', cmap='Blues', ax=ax)
            ax.set_xlabel('Gender')
            ax.set_ylabel('Virus Type')
            st.pyplot(fig)
        else:
            st.warning("Gender column ('Cinsiyet') not found or only one gender present for virus distribution.")

if "Time Series Analysis" in analysis_options and 'Analyzed_Date' in df.columns and actual_virus_columns and df['Analyzed_Date'].notna().any():
    st.header("3. Time Series Analysis")

    st.subheader("Virus Trends Over Time")
    time_option = st.selectbox("Time granularity", ["Monthly", "Weekly", "Yearly"])

    if time_option == "Monthly":
        monthly_virus = df.groupby(['Year', 'Month'])[actual_virus_columns].sum().reset_index()
        monthly_virus['Year_Month'] = pd.to_datetime(
            monthly_virus['Year'].astype(str) + '-' + monthly_virus['Month'].astype(str) + '-01'
        )
        monthly_virus = monthly_virus.sort_values('Year_Month')

        if not monthly_virus.empty:
            fig = px.line(monthly_virus, x='Year_Month', y=actual_virus_columns,
                          title='Monthly Virus Trends')
            st.plotly_chart(fig, use_container_width=True)

            fig = px.line(monthly_virus, x='Month', y=actual_virus_columns,
                          facet_col='Year', facet_col_wrap=4,
                          title='Monthly Virus Trends by Year')
            fig.update_xaxes(nticks=12)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No monthly virus data available for plotting.")

    elif time_option == "Weekly":
        if 'Week_of_Year' in df.columns:
            weekly_virus = df.groupby(['Year', 'Week_of_Year'])[actual_virus_columns].sum().reset_index()
            if not weekly_virus.empty:
                fig = px.line(weekly_virus, x='Week_of_Year', y=actual_virus_columns,
                              facet_col='Year', facet_col_wrap=4,
                              title='Weekly Virus Trends by Year')
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.warning("No weekly virus data available for plotting.")
        else:
            st.warning("Weekly data not available. 'Week_of_Year' column is missing or could not be generated.")

    elif time_option == "Yearly":
        yearly_virus = df.groupby('Year')[actual_virus_columns].sum().reset_index()
        if not yearly_virus.empty:
            fig = px.bar(yearly_virus, x='Year', y=actual_virus_columns,
                          title='Yearly Virus Trends')
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.warning("No yearly virus data available for plotting.")
else:
    if "Time Series Analysis" in analysis_options:
        st.warning("Time series analysis skipped: 'Analyzed_Date' column is missing or contains no valid dates, or no virus columns are identified.")


if "Environmental Factors" in analysis_options and env_temp_col in df.columns and env_hum_col in df.columns and actual_virus_columns:
    # Ensure there's enough non-NA data for plotting
    if df[env_temp_col].notna().any() and df[env_hum_col].notna().any():
        st.header("4. Environmental Factors Analysis")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Temperature and Humidity Distribution")
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
            sns.histplot(df[env_temp_col].dropna(), bins=20, kde=True, color='crimson', ax=ax1) # dropna for histplot
            ax1.set_title('Monthly Temperature')
            ax1.set_xlabel('Temperature (°C)')

            sns.histplot(df[env_hum_col].dropna(), bins=20, kde=True, color='royalblue', ax=ax2) # dropna for histplot
            ax2.set_title('Monthly Humidity')
            ax2.set_xlabel('Humidity (%)')

            plt.tight_layout()
            st.pyplot(fig)

        with col2:
            st.subheader("Virus Distribution by Temperature")
            # Ensure enough bins are present after cutting
            temp_bins = pd.cut(df[env_temp_col], bins=6)
            if not temp_bins.empty and temp_bins.nunique() > 1: # Check for multiple bins
                temp_virus = df.groupby(temp_bins)[actual_virus_columns].sum().T
                if not temp_virus.empty:
                    fig, ax = plt.subplots(figsize=(10, 8))
                    sns.heatmap(temp_virus, annot=True, fmt='d', cmap='coolwarm', ax=ax)
                    ax.set_xlabel('Temperature Range (°C)')
                    ax.set_ylabel('Virus Type')
                    st.pyplot(fig)
                else:
                    st.warning("Not enough data to group viruses by temperature ranges.")
            else:
                st.warning("No valid temperature data to create ranges or only one range present.")

        st.subheader("Virus Distribution by Humidity")
        humidity_bins = pd.cut(df[env_hum_col], bins=6)
        if not humidity_bins.empty and humidity_bins.nunique() > 1: # Check for multiple bins
            humidity_virus = df.groupby(humidity_bins)[actual_virus_columns].sum().T
            if not humidity_virus.empty:
                fig, ax = plt.subplots(figsize=(10, 8))
                sns.heatmap(humidity_virus, annot=True, fmt='d', cmap='summer', ax=ax)
                ax.set_xlabel('Humidity Range (%)')
                ax.set_ylabel('Virus Type')
                st.pyplot(fig)
            else:
                st.warning("Not enough data to group viruses by humidity ranges.")
        else:
            st.warning("No valid humidity data to create ranges or only one range present.")
    else:
        st.warning("Environmental Factors Analysis skipped: Temperature or Humidity columns are missing or contain no valid data.")
else:
    if "Environmental Factors" in analysis_options:
        st.warning("Environmental Factors Analysis skipped: Temperature, Humidity, or virus columns are not found.")

if "Co-infections" in analysis_options and actual_virus_columns and 'Total_Virus_Count' in df.columns:
    st.header("5. Co-infections Analysis")
    co_infected = df[df['Total_Virus_Count'] > 1]
    st.metric("Patients with co-infections",
              f"{len(co_infected)} ({(len(co_infected)/len(df)*100):.2f}%)")

    try:
        from mlxtend.frequent_patterns import apriori, association_rules

        st.subheader("Virus Association Rules")

        virus_binary = df[actual_virus_columns].applymap(lambda x: 1 if x > 0 else 0)
        virus_binary = virus_binary[virus_binary.sum(axis=1) > 0] # Filter out rows with no viruses

        if not virus_binary.empty:
            # Check if 'min_support' is defined, otherwise use a default
            min_support_val = min_support if 'min_support' in locals() else 0.005
            frequent_itemsets = apriori(virus_binary, min_support=min_support_val, use_colnames=True)

            if not frequent_itemsets.empty:
                rules = association_rules(frequent_itemsets, metric="lift", min_threshold=1)
                st.dataframe(rules.sort_values('lift', ascending=False).head(10))
            else:
                st.warning("No frequent itemsets found with current support threshold. Try lowering the 'Minimum support' in Advanced Options.")
        else:
            st.warning("No valid data for association rule mining. Ensure there are patients with at least one virus detected.")
    except ImportError:
        st.error("mlxtend library not installed. Please install with: pip install mlxtend")
    except Exception as e:
        st.error(f"An error occurred during association rule mining: {e}")

# --- MODIFIED SECTION: Advanced Statistical Analysis (formerly 9, now replaces 6) ---
if "Advanced Statistical Analysis" in analysis_options and actual_virus_columns:
    st.header("6. Advanced Statistical Analysis (Interactive)")
    st.write("Explore statistical relationships between a selected virus and other demographic/environmental factors.")

    col_stat_1, col_stat_2 = st.columns(2)

    with col_stat_1:
        selected_virus_stat = st.selectbox(
            "Select a Virus to Analyze:",
            actual_virus_columns,
            key='selected_virus_stat'
        )

    with col_stat_2:
        # Define available grouping variables based on what's present in the DataFrame
        available_grouping_vars = []
        if 'Cinsiyet' in df.columns and df['Cinsiyet'].nunique(dropna=True) > 1:
            available_grouping_vars.append('Cinsiyet')
        if 'Yaş_Grubu' in df.columns and df['Yaş_Grubu'].nunique(dropna=True) > 1:
            available_grouping_vars.append('Yaş_Grubu')
        if 'Yaş' in df.columns and df['Yaş'].notna().any():
            available_grouping_vars.append('Yaş')
        if env_temp_col in df.columns and df[env_temp_col].notna().any():
            available_grouping_vars.append(env_temp_col)
        if env_hum_col in df.columns and df[env_hum_col].notna().any():
            available_grouping_vars.append(env_hum_col)
        if 'Month' in df.columns and df['Month'].nunique(dropna=True) > 1:
            available_grouping_vars.append('Month')
        if 'Season' in df.columns and df['Season'].nunique(dropna=True) > 1:
            available_grouping_vars.append('Season')

        if not available_grouping_vars:
            st.warning("No suitable grouping variables found for advanced statistical analysis (requires 'Cinsiyet', 'Yaş_Grubu', 'Yaş', 'Temp_ay_ortalaması', 'Humidity_ayortalaması', 'Month', or 'Season' with sufficient data).")
        else:
            selected_grouping_var = st.selectbox(
                "Select a Grouping Variable:",
                available_grouping_vars,
                key='selected_grouping_var'
            )

    if selected_virus_stat and selected_grouping_var:
        st.subheader(f"Analysis of {selected_virus_stat} by {selected_grouping_var}")

        # Create a binary target for the selected virus (0 for negative, 1 for positive)
        df_analysis = df.copy()
        df_analysis[f'{selected_virus_stat}_Present'] = (df_analysis[selected_virus_stat] > 0).astype(int)

        # Drop rows where either the virus presence or the grouping variable is NaN
        df_analysis.dropna(subset=[f'{selected_virus_stat}_Present', selected_grouping_var], inplace=True)

        if df_analysis.empty:
            st.warning("No data available for the selected virus and grouping variable after filtering missing values.")
        elif df_analysis[f'{selected_virus_stat}_Present'].nunique() < 2:
            st.warning(f"The selected virus '{selected_virus_stat}' has only one outcome (all present or all absent) after filtering. Statistical tests cannot be performed.")
        else:
            # Determine if the grouping variable is categorical or continuous
            is_categorical = selected_grouping_var in ['Cinsiyet', 'Yaş_Grubu', 'Month', 'Season']

            if is_categorical:
                # Chi-square test for categorical grouping variables
                st.write("Performing Chi-square test for independence between a categorical variable and virus presence.")
                contingency_table = pd.crosstab(df_analysis[selected_grouping_var], df_analysis[f'{selected_virus_stat}_Present'])

                # Ensure sufficient data for chi-square (at least 2x2 and no zero expected frequencies)
                if contingency_table.shape[0] > 1 and contingency_table.shape[1] > 1 and (contingency_table.min().min() > 0 or (contingency_table.sum(axis=1) > 0).all() and (contingency_table.sum(axis=0) > 0).all()):
                    try:
                        chi2, p, _, expected = stats.chi2_contingency(contingency_table)
                        st.metric(f"Chi-square test p-value for {selected_virus_stat} vs {selected_grouping_var}", f"{p:.4f}")

                        if p < 0.05:
                            st.success(f"**Interpretation:** There is a statistically significant association between **{selected_grouping_var}** and **{selected_virus_stat}** presence (p < 0.05). This suggests that the distribution of {selected_virus_stat} presence is not independent of {selected_grouping_var}. For example, certain {selected_grouping_var} categories might have a higher or lower proportion of {selected_virus_stat} cases.")
                        else:
                            st.info(f"**Interpretation:** No statistically significant association found between **{selected_grouping_var}** and **{selected_virus_stat}** presence (p >= 0.05). This suggests that the distribution of {selected_virus_stat} presence is independent of {selected_grouping_var}. The proportion of {selected_virus_stat} cases is similar across different {selected_grouping_var} categories.")

                        st.write("Contingency Table:")
                        st.dataframe(contingency_table)

                        fig, ax = plt.subplots(figsize=(10, 6))
                        sns.countplot(data=df_analysis, x=selected_grouping_var, hue=f'{selected_virus_stat}_Present', palette='viridis', ax=ax)
                        ax.set_title(f'Distribution of {selected_virus_stat} by {selected_grouping_var}')
                        ax.set_xlabel(selected_grouping_var)
                        ax.set_ylabel('Count')
                        st.pyplot(fig)
                    except ValueError as ve:
                        st.warning(f"Could not perform Chi-square test: {ve}. This often happens if expected frequencies are too low in some cells. Consider aggregating categories if possible.")
                else:
                    st.warning(f"Insufficient data or imbalanced categories for Chi-square test for {selected_virus_stat} vs {selected_grouping_var}. Ensure there are at least two unique values in both the virus presence and grouping variable, and sufficient observations in each category.")

            else: # Continuous grouping variable (Yaş, Temp_ay_ortalaması, Humidity_ayortalaması)
                # Perform Independent t-test (comparing means of continuous variable for virus present vs. absent)
                st.write("Performing Independent t-test to compare the mean of a continuous variable between virus present and absent groups.")
                group_present = df_analysis[df_analysis[f'{selected_virus_stat}_Present'] == 1][selected_grouping_var].dropna()
                group_absent = df_analysis[df_analysis[f'{selected_virus_stat}_Present'] == 0][selected_grouping_var].dropna()

                if len(group_present) > 1 and len(group_absent) > 1 and group_present.var() > 0 and group_absent.var() > 0:
                    t_stat, p_value = stats.ttest_ind(group_present, group_absent, equal_var=False) # Welch's t-test
                    st.metric(f"Independent t-test p-value for {selected_virus_stat} vs {selected_grouping_var}", f"{p_value:.4g}")
                    if p_value < 0.05:
                        st.success(f"**Interpretation:** There is a statistically significant difference in the mean **{selected_grouping_var}** between patients with and without **{selected_virus_stat}** (p < 0.05). This suggests that the average {selected_grouping_var} value is different for individuals who test positive for {selected_virus_stat} compared to those who test negative.")
                    else:
                        st.info(f"**Interpretation:** No statistically significant difference found in the mean **{selected_grouping_var}** between patients with and without **{selected_virus_stat}** (p >= 0.05). This suggests that the average {selected_grouping_var} value is similar for both groups.")

                    st.write(f"Mean {selected_grouping_var} for {selected_virus_stat} Present: **{group_present.mean():.2f}**")
                    st.write(f"Mean {selected_grouping_var} for {selected_virus_stat} Absent: **{group_absent.mean():.2f}**")

                    fig, ax = plt.subplots(figsize=(10, 6))
                    sns.boxplot(data=df_analysis, x=f'{selected_virus_stat}_Present', y=selected_grouping_var, palette='pastel', ax=ax)
                    ax.set_title(f'{selected_grouping_var} Distribution by {selected_virus_stat} Presence')
                    ax.set_xlabel(f'{selected_virus_stat} Present (0=No, 1=Yes)')
                    ax.set_ylabel(selected_grouping_var)
                    st.pyplot(fig)

                    fig_kde, ax_kde = plt.subplots(figsize=(10, 6))
                    sns.kdeplot(data=df_analysis, x=selected_grouping_var, hue=f'{selected_virus_stat}_Present', fill=True, common_norm=False, palette='coolwarm', ax=ax_kde)
                    ax_kde.set_title(f'Density Plot of {selected_grouping_var} by {selected_virus_stat} Presence')
                    ax_kde.set_xlabel(selected_grouping_var)
                    ax_kde.set_ylabel('Density')
                    st.pyplot(fig_kde)

                    # Also add correlation for continuous variables
                    correlation = df_analysis[[selected_grouping_var, f'{selected_virus_stat}_Present']].corr().iloc[0, 1]
                    st.metric(f"Pearson Correlation between {selected_grouping_var} and {selected_virus_stat} presence", f"{correlation:.4f}")
                    if abs(correlation) > 0.1: # A simple threshold for 'notable' correlation
                        st.info(f"**Interpretation:** There is a notable correlation of {correlation:.2f} between {selected_grouping_var} and {selected_virus_stat} presence. A positive correlation means as {selected_grouping_var} increases, {selected_virus_stat} presence tends to increase (and vice-versa for negative correlation).")
                    else:
                        st.info(f"**Interpretation:** The correlation between {selected_grouping_var} and {selected_virus_stat} presence is weak ({correlation:.2f}). This suggests a limited linear relationship between these two variables.")

                else:
                    st.warning(f"Insufficient data or no variance in groups for t-test for {selected_virus_stat} vs {selected_grouping_var}.")
    else:
        st.info("Please select a virus and a grouping variable to perform advanced statistical analysis.")

st.header("7. Logistic Regression Analysis")
st.write("""
Çok değişkenli (multivariable) lojistik regresyon ile seçtiğiniz virüsün pozitifliği üzerinde
yaş, cinsiyet, sıcaklık, nem, sezon vb. değişkenlerin etkisini inceleyebilirsiniz.
""")

# Hedef virüs seçimi
target_virus = st.selectbox(
    "Hedef (Dependent Variable): Analiz etmek istediğiniz virüsü seçin.",
    actual_virus_columns,
    key="logreg_target_virus"
)

# Bağımsız değişkenleri belirleyelim
potential_predictors = []
if 'Yaş' in df.columns: potential_predictors.append('Yaş')
if 'Cinsiyet' in df.columns: potential_predictors.append('Cinsiyet')
if env_temp_col in df.columns: potential_predictors.append(env_temp_col)
if env_hum_col in df.columns: potential_predictors.append(env_hum_col)
if 'Season' in df.columns: potential_predictors.append('Season')

selected_predictors = st.multiselect(
    "Bağımsız değişkenler (Predictors): İstediğiniz kadar seçebilirsiniz.",
    potential_predictors,
    default=potential_predictors
)

if st.button("Lojistik Regresyonu Çalıştır"):
    df_logreg = df.copy()
    # Hedefi binary'ye çevir
    df_logreg["target"] = (df_logreg[target_virus] > 0).astype(int)
    used_cols = ['target'] + selected_predictors
    df_logreg = df_logreg[used_cols].dropna()
    
    # Kategorik değişkenleri dummy'ye çevir
    categorical_vars = [col for col in selected_predictors if
                        str(df_logreg[col].dtype) == 'object' or col in ['Cinsiyet', 'Season']]
    df_logreg = pd.get_dummies(df_logreg, columns=categorical_vars, drop_first=True)
    df_logreg = df_logreg.astype(float)

    if df_logreg.shape[0] < 10:
        st.warning("Regresyon için yeterli veri yok!")
    elif df_logreg['target'].nunique() < 2:
        st.warning("Seçtiğiniz virüs için yeterli pozitif/negatif vaka yok!")
    else:
        X = df_logreg.drop(columns=['target']).astype(float)
        y = df_logreg['target'].astype(float)
        X = sm.add_constant(X)

        logit_model = sm.Logit(y, X)
        try:
            result = logit_model.fit(disp=0)
        except Exception as e:
            st.error(f"Regresyon sırasında hata: {e}")
            result = None

        if result:
            st.subheader("Lojistik Regresyon Sonuçları")
            summary_df = result.summary2().tables[1]
            summary_df["Odds Ratio"] = summary_df["Coef."].apply(np.exp)
            summary_df["CI Lower"] = (summary_df["Coef."] - 1.96*summary_df["Std.Err."]).apply(np.exp)
            summary_df["CI Upper"] = (summary_df["Coef."] + 1.96*summary_df["Std.Err."]).apply(np.exp)

            # p<0.05 olanları vurgulayalım
            def highlight_sig(row):
                # row burada bir pandas Series, yani bir satır
                p = row["P>|z|"]
                if p < 0.05:
                    return ['font-weight: bold; background-color: #d0f5dd'] * len(row)
                else:
                    return [''] * len(row)
            
            st.dataframe(
                summary_df[['Odds Ratio', 'CI Lower', 'CI Upper', 'P>|z|']],
                use_container_width=True
            )

            # === Odds Ratio Forest Plot ===
            summary_plot = summary_df.drop("const", errors="ignore")
            fig, ax = plt.subplots(figsize=(7, 0.7 * len(summary_plot)))
            ax.errorbar(
                summary_plot["Odds Ratio"],
                summary_plot.index,
                xerr=[
                    summary_plot["Odds Ratio"] - summary_plot["CI Lower"],
                    summary_plot["CI Upper"] - summary_plot["Odds Ratio"]
                ],
                fmt='o', color='teal', capsize=5, label='Odds Ratio (95% CI)'
            )
            ax.axvline(1, color='grey', linestyle='--', lw=1)
            ax.set_xlabel("Odds Ratio (OR, log scale)")
            ax.set_xscale("log")
            ax.set_title("Odds Ratio ve 95% Güven Aralığı (Lojistik Regresyon)")
            ax.set_ylabel("Değişkenler")
            plt.tight_layout()
            st.pyplot(fig)

            # === ROC Curve ===
            y_pred_prob = result.predict(X)
            fpr, tpr, _ = roc_curve(y, y_pred_prob)
            auc_score = roc_auc_score(y, y_pred_prob)
            fig2, ax2 = plt.subplots()
            ax2.plot(fpr, tpr, label=f'AUC = {auc_score:.2f}')
            ax2.plot([0, 1], [0, 1], '--', color='grey')
            ax2.set_xlabel("1 - Özgüllük (False Positive Rate)")
            ax2.set_ylabel("Duyarlılık (True Positive Rate)")
            ax2.set_title("ROC Curve")
            ax2.legend()
            st.pyplot(fig2)

            # === Modelin Genel p-Değeri (Likelihood Ratio Test) ===
            llf_full = result.llf
            llf_null = result.llnull
            df_full = result.df_model + 1  # +1: sabit/intercept dahil
            df_null = 1  # Sadece sabitli model
            lr_stat = 2 * (llf_full - llf_null)
            lr_df = df_full - df_null
            lr_pvalue = sps.chi2.sf(lr_stat, lr_df)

            st.markdown(f"""
                **Modelin genel anlamlılığı (Likelihood Ratio test):**  
                - Test statistic: `{lr_stat:.2f}`
                - df: `{lr_df}`
                - p-value: `{lr_pvalue:.4g}`
                """)


            if lr_pvalue < 0.05:
                st.success("Model genel olarak anlamlıdır (p < 0.05). Yani bağımsız değişkenlerin tamamı bir arada, hedef değişkenin açıklanmasında anlamlı katkı sağlıyor.")
            else:
                st.warning("Modelin genel anlamlılığı yetersiz (p ≥ 0.05). Modelin açıklama gücü zayıf olabilir.")

            st.markdown("""
            **Odds Ratio (OR) > 1:** Değişken arttıkça hedef virüs pozitifliği olasılığı artıyor.<br>
            **OR < 1:** Değişken arttıkça olasılık azalıyor.<br>
            **p < 0.05:** İstatistiksel olarak anlamlı.<br>
            <br>
            *Sonuçlar model varsayımlarına ve verinin yapısına göre yorumlanmalıdır.*
            """, unsafe_allow_html=True)


if "Clustering" in analysis_options:
    st.header("7. Clustering Analysis") # Re-numbered to 7
    st.write("Group patients based on selected features to identify distinct patterns.")

    # Get all numerical columns for clustering
    numerical_cols = df.select_dtypes(include=np.number).columns.tolist()
    # Exclude 'Total_Virus_Count', 'Year', 'Month', 'Week_of_Year' as they are derived or time-based
    # and might not be ideal for direct clustering of patient characteristics.
    # Include 'Yaş', env_temp_col, env_hum_col and all actual virus columns as primary candidates.
    clustering_feature_candidates = ['Yaş', env_temp_col, env_hum_col] + [col for col in actual_virus_columns if col in numerical_cols]
    clustering_feature_candidates = list(set(clustering_feature_candidates)) # Remove duplicates

    if not clustering_feature_candidates:
        st.warning("No suitable numerical columns found for clustering. Please ensure 'Yaş', temperature, humidity, or virus columns are present and numerical.")
    else:
        st.subheader("Select Features for Clustering")
        col_cluster_feat_1, col_cluster_feat_2 = st.columns(2)
        with col_cluster_feat_1:
            selected_x_feature = st.selectbox(
                "Select X-axis feature for visualization:",
                clustering_feature_candidates,
                index=0 if 'Yaş' in clustering_feature_candidates else (clustering_feature_candidates.index(env_temp_col) if env_temp_col in clustering_feature_candidates else 0), # Default to 'Yaş' or first available
                key="cluster_x_feature"
            )
        with col_cluster_feat_2:
            # Ensure y-feature is different from x-feature
            y_options = [f for f in clustering_feature_candidates if f != selected_x_feature]
            selected_y_feature = st.selectbox(
                "Select Y-axis feature for visualization:",
                y_options,
                index=0 if env_temp_col in y_options else (y_options.index('Yaş') if 'Yaş' in y_options else 0), # Default to Temp or first available different from x
                key="cluster_y_feature"
            )

        # Prepare data for clustering based on selected features
        # This data is used only for the 2D scatter plot visualization
        cluster_data_for_plot = df[[selected_x_feature, selected_y_feature]].copy()
        cluster_data_for_plot[selected_x_feature] = pd.to_numeric(cluster_data_for_plot[selected_x_feature], errors='coerce')
        cluster_data_for_plot[selected_y_feature] = pd.to_numeric(cluster_data_for_plot[selected_y_feature], errors='coerce')
        cluster_data_for_plot.dropna(inplace=True)

        # Use all relevant numerical features for the actual clustering algorithm
        # This ensures the algorithm uses more information than just the two selected for visualization
        all_numerical_for_clustering = [col for col in numerical_cols if col not in ['Total_Virus_Count', 'Year', 'Month', 'Week_of_Year']]
        cluster_data_for_algorithm = df[all_numerical_for_clustering].copy()
        for col in all_numerical_for_clustering:
            cluster_data_for_algorithm[col] = pd.to_numeric(cluster_data_for_algorithm[col], errors='coerce')
        cluster_data_for_algorithm.dropna(inplace=True)


        if cluster_data_for_algorithm.empty or cluster_data_for_plot.empty:
            st.warning("No valid data for clustering after dropping missing values in selected columns or all numerical columns.")
        else:
            clustering_method = st.selectbox(
                "Select Clustering Method:",
                ["K-Means", "DBSCAN"],
                key="clustering_method_selection"
            )

            df_copy_for_cluster = df.copy() # Create a copy to add cluster labels

            if clustering_method == "K-Means":
                st.subheader("K-Means Clustering")
                st.write("K-Means aims to partition n observations into k clusters in which each observation belongs to the cluster with the nearest mean (cluster centroids).")
                st.write("It works best with spherical clusters of similar size and density.")

                num_clusters = st.number_input("Number of clusters (k)", 2, 10, 4, key="kmeans_n_clusters")
                kmeans_random_state = st.number_input("Random state for K-Means", 42, key="kmeans_random_state")

                if len(cluster_data_for_algorithm) >= num_clusters:
                    kmeans = KMeans(
                        n_clusters=num_clusters,
                        random_state=kmeans_random_state,
                        n_init='auto'
                    )
                    # Apply clustering on the full set of numerical features
                    df_copy_for_cluster.loc[cluster_data_for_algorithm.index, 'Cluster'] = kmeans.fit_predict(cluster_data_for_algorithm)
                    df_copy_for_cluster['Cluster'] = df_copy_for_cluster['Cluster'].astype(str) # For plotting

                    st.subheader(f"Cluster Visualization ({selected_x_feature} vs {selected_y_feature})")
                    # Plot using the selected two features for visualization
                    fig = px.scatter(
                        df_copy_for_cluster.dropna(subset=['Cluster', selected_x_feature, selected_y_feature]),
                        x=selected_x_feature, y=selected_y_feature, color='Cluster',
                        title=f'Patient Clusters ({selected_x_feature} vs {selected_y_feature})',
                        labels={selected_x_feature: selected_x_feature, selected_y_feature: selected_y_feature}
                    )
                    st.plotly_chart(fig, use_container_width=True)

                    if 'Cluster' in df_copy_for_cluster.columns and actual_virus_columns:
                        st.subheader("Average Virus Presence by Cluster")
                        # Calculate mean presence (0 or 1) for each virus in each cluster
                        virus_profile_by_cluster = df_copy_for_cluster.groupby('Cluster')[actual_virus_columns].mean().T
                        if not virus_profile_by_cluster.empty:
                            fig, ax = plt.subplots(figsize=(12, 8))
                            sns.heatmap(virus_profile_by_cluster, annot=True, fmt=".2f", cmap='Greens', ax=ax)
                            ax.set_title('Average Virus Presence in Each Cluster')
                            ax.set_xlabel('Cluster')
                            ax.set_ylabel('Virus Type')
                            st.pyplot(fig)
                        else:
                            st.warning("No virus data available for clusters to show average presence.")

                    st.subheader("Cluster Characteristics by Other Features")
                    # Display distributions of other relevant features for each cluster
                    characteristics_features = ['Yaş', 'Cinsiyet', env_temp_col, env_hum_col, 'Month', 'Season']
                    available_characteristics = [f for f in characteristics_features if f in df_copy_for_cluster.columns and df_copy_for_cluster[f].notna().any()]

                    for feature in available_characteristics:
                        if feature in ['Yaş', env_temp_col, env_hum_col]: # Continuous
                            fig, ax = plt.subplots(figsize=(10, 6))
                            sns.boxplot(data=df_copy_for_cluster, x='Cluster', y=feature, palette='pastel', ax=ax)
                            ax.set_title(f'{feature} Distribution by Cluster')
                            st.pyplot(fig)
                        elif feature in ['Cinsiyet', 'Month', 'Season']: # Categorical
                            fig, ax = plt.subplots(figsize=(10, 6))
                            sns.countplot(data=df_copy_for_cluster, x='Cluster', hue=feature, palette='viridis', ax=ax)
                            ax.set_title(f'{feature} Distribution by Cluster')
                            ax.tick_params(axis='x', rotation=45)
                            st.pyplot(fig)
                else:
                    st.warning(f"Insufficient data for K-Means clustering (less than {num_clusters} samples). Please check your data or reduce the number of clusters.")

            elif clustering_method == "DBSCAN":
                st.subheader("DBSCAN Clustering")
                st.write("DBSCAN (Density-Based Spatial Clustering of Applications with Noise) groups together points that are closely packed together, marking as outliers points that lie alone in low-density regions. It's good for finding clusters of arbitrary shape and identifying noise.")
                st.write("**Parameters:**")
                st.markdown(
                    """
                    - **Epsilon (eps):** The maximum distance between two samples for one to be considered as in the neighborhood of the other.
                    - **Min Samples (min_samples):** The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.
                    """
                )

                eps_val = st.slider("Epsilon (eps)", 0.1, 10.0, 0.5, 0.1, key="dbscan_eps")
                min_samples_val = st.number_input("Min Samples (min_samples)", 1, 50, 5, key="dbscan_min_samples")

                if len(cluster_data_for_algorithm) > 0:
                    dbscan = DBSCAN(eps=eps_val, min_samples=min_samples_val)
                    # Apply clustering on the full set of numerical features
                    dbscan_labels = dbscan.fit_predict(cluster_data_for_algorithm)
                    df_copy_for_cluster.loc[cluster_data_for_algorithm.index, 'Cluster'] = dbscan_labels
                    df_copy_for_cluster['Cluster'] = df_copy_for_cluster['Cluster'].astype(str) # For plotting

                    n_clusters_ = len(set(dbscan_labels)) - (1 if -1 in dbscan_labels else 0)
                    n_noise_ = list(dbscan_labels).count(-1)

                    st.info(f"Estimated number of clusters: {n_clusters_}")
                    st.info(f"Estimated number of noise points: {n_noise_}")

                    st.subheader(f"Cluster Visualization ({selected_x_feature} vs {selected_y_feature})")
                    # Plot using the selected two features for visualization
                    fig = px.scatter(
                        df_copy_for_cluster.dropna(subset=['Cluster', selected_x_feature, selected_y_feature]),
                        x=selected_x_feature, y=selected_y_feature, color='Cluster',
                        title=f'Patient Clusters ({selected_x_feature} vs {selected_y_feature}) - DBSCAN',
                        labels={selected_x_feature: selected_x_feature, selected_y_feature: selected_y_feature},
                        color_discrete_map={'-1': 'black'} # Map noise to black
                    )
                    st.plotly_chart(fig, use_container_width=True)

                    if 'Cluster' in df_copy_for_cluster.columns and actual_virus_columns:
                        st.subheader("Average Virus Presence by Cluster (Excluding Noise)")
                        # Filter out noise points (-1) for this heatmap
                        cluster_virus = df_copy_for_cluster[df_copy_for_cluster['Cluster'] != '-1'].groupby('Cluster')[actual_virus_columns].mean().T
                        if not cluster_virus.empty:
                            fig, ax = plt.subplots(figsize=(12, 8))
                            sns.heatmap(cluster_virus, annot=True, fmt=".2f", cmap='Greens', ax=ax)
                            ax.set_title('Average Virus Presence in Each Cluster (Excluding Noise)')
                            ax.set_xlabel('Cluster')
                            ax.set_ylabel('Virus Type')
                            st.pyplot(fig)
                        else:
                            st.warning("No virus data available for non-noise clusters to show average presence.")

                    st.subheader("Cluster Characteristics by Other Features (Excluding Noise)")
                    characteristics_features = ['Yaş', 'Cinsiyet', env_temp_col, env_hum_col, 'Month', 'Season']
                    # Filter out noise points for these plots
                    df_filtered_for_dbscan_plots = df_copy_for_cluster[df_copy_for_cluster['Cluster'] != '-1'].copy()
                    available_characteristics = [f for f in characteristics_features if f in df_filtered_for_dbscan_plots.columns and df_filtered_for_dbscan_plots[f].notna().any()]

                    for feature in available_characteristics:
                        if feature in ['Yaş', env_temp_col, env_hum_col]: # Continuous
                            fig, ax = plt.subplots(figsize=(10, 6))
                            sns.boxplot(data=df_filtered_for_dbscan_plots, x='Cluster', y=feature, palette='pastel', ax=ax)
                            ax.set_title(f'{feature} Distribution by Cluster (Excluding Noise)')
                            st.pyplot(fig)
                        elif feature in ['Cinsiyet', 'Month', 'Season']: # Categorical
                            fig, ax = plt.subplots(figsize=(10, 6))
                            sns.countplot(data=df_filtered_for_dbscan_plots, x='Cluster', hue=feature, palette='viridis', ax=ax)
                            ax.set_title(f'{feature} Distribution by Cluster (Excluding Noise)')
                            ax.tick_params(axis='x', rotation=45)
                            st.pyplot(fig)
                else:
                    st.warning("Insufficient data for DBSCAN clustering.")
            else:
              st.warning("Clustering analysis skipped due to missing required numerical columns or empty data after preprocessing.")


# Function to generate HTML report (adapted from previous version, make sure it's defined)
def generate_html_report_streamlit(dataframe, temp_directory, virus_cols, age_labels, env_temp_col_name, env_hum_col_name):
    """Generates an HTML report for the analysis."""
    report_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Respiratory Virus Analysis Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; line-height: 1.6; margin: 20px; }}
            h1, h2, h3 {{ color: #2c3e50; border-bottom: 1px solid #ddd; padding-bottom: 5px; }}
            img {{ max-width: 100%; height: auto; display: block; margin: 10px 0; }}
            table {{ width: 100%; border-collapse: collapse; margin-bottom: 20px; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <h1>Respiratory Virus Analysis Report</h1>
        <p>Report generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>Total records analyzed: {len(dataframe)}</p>

        <h2>1. Demographic Analysis</h2>
    """
    # Define static graph directory inside the temporary directory
    static_graph_dir = os.path.join(temp_directory, "static_graphs")
    os.makedirs(static_graph_dir, exist_ok=True) # Ensure it exists

    # Save all figures for the report
    # Demographic Analysis
    if 'Cinsiyet' in dataframe.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        dataframe['Cinsiyet'].value_counts().plot(kind='pie', autopct='%1.1f%%', colors=['skyblue', 'lightcoral'], ax=ax)
        ax.set_ylabel('')
        fig.savefig(os.path.join(static_graph_dir, "demografi_cinsiyet.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        report_content += """
        <h3>Gender Distribution</h3>
        <img src="static_graphs/demografi_cinsiyet.png" alt="Gender Distribution">
        """

    if 'Yaş' in dataframe.columns:
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.histplot(dataframe['Yaş'], bins=30, kde=True, color='teal', ax=ax)
        ax.set_title('Age Distribution')
        ax.set_xlabel('Age')
        ax.set_ylabel('Patient Count')
        fig.savefig(os.path.join(static_graph_dir, "yas_dagilimi.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        report_content += """
        <h3>Age Distribution</h3>
        <img src="static_graphs/yas_dagilimi.png" alt="Age Distribution">
        """

        if 'Yaş_Grubu' in dataframe.columns:
            fig, ax = plt.subplots(figsize=(10, 6))
            dataframe['Yaş_Grubu'].value_counts().sort_index().plot(kind='bar', color='mediumseagreen', ax=ax)
            ax.set_title('Age Group Distribution')
            ax.set_xlabel('Age Group')
            ax.set_ylabel('Patient Count')
            ax.tick_params(axis='x', rotation=45)
            fig.savefig(os.path.join(static_graph_dir, "yas_grubu_dagilimi.png"), bbox_inches='tight', dpi=300)
            plt.close(fig)
            report_content += """
            <h3>Age Group Distribution</h3>
            <img src="static_graphs/yas_grubu_dagilimi.png" alt="Age Group Distribution">
            """

    # Virus Distribution
    if virus_cols:
        report_content += """
        <h2>2. Virus Distribution Analysis</h2>
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        dataframe[virus_cols].sum().sort_values(ascending=False).plot(kind='barh', color='darkorange', ax=ax)
        ax.set_title('Virus Case Counts')
        ax.set_xlabel('Case Count')
        ax.set_ylabel('Virus Type')
        fig.savefig(os.path.join(static_graph_dir, "virus_dagilimi.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        report_content += """
        <h3>Virus Case Counts</h3>
        <img src="static_graphs/virus_dagilimi.png" alt="Virus Case Counts">
        """

        if 'Yaş_Grubu' in dataframe.columns:
            fig, ax = plt.subplots(figsize=(10, 8))
            age_virus = dataframe.groupby('Yaş_Grubu')[virus_cols].sum().T
            sns.heatmap(age_virus, annot=True, fmt='d', cmap='YlOrRd', ax=ax)
            ax.set_title('Virus Distribution by Age Group')
            ax.set_xlabel('Age Group')
            ax.set_ylabel('Virus Type')
            fig.savefig(os.path.join(static_graph_dir, "yas_grubu_virus.png"), bbox_inches='tight', dpi=300)
            plt.close(fig)
            report_content += """
            <h3>Virus Distribution by Age Group</h3>
            <img src="static_graphs/yas_grubu_virus.png" alt="Virus Distribution by Age Group">
            """

        if 'Cinsiyet' in dataframe.columns:
            fig, ax = plt.subplots(figsize=(10, 8))
            gender_virus = dataframe.groupby('Cinsiyet')[virus_cols].sum().T
            sns.heatmap(gender_virus, annot=True, fmt='d', cmap='Blues', ax=ax)
            ax.set_title('Virus Distribution by Gender')
            ax.set_xlabel('Gender')
            ax.set_ylabel('Virus Type')
            fig.savefig(os.path.join(static_graph_dir, "cinsiyet_virus.png"), bbox_inches='tight', dpi=300)
            plt.close(fig)
            report_content += """
            <h3>Virus Distribution by Gender</h3>
            <img src="static_graphs/cinsiyet_virus.png" alt="Virus Distribution by Gender">
            """

    # Time Series Analysis
    if 'Analyzed_Date' in dataframe.columns and virus_cols and dataframe['Analyzed_Date'].notna().any():
        report_content += """
        <h2>3. Time Series Analysis</h2>
        """
        fig, ax = plt.subplots(figsize=(14, 8))
        month_virus = dataframe.groupby(dataframe['Analyzed_Date'].dt.month)[virus_cols].sum().T
        sns.heatmap(month_virus, annot=True, fmt='d', cmap='viridis', ax=ax)
        ax.set_title('Monthly Virus Distribution (All Years Total)')
        ax.set_xlabel('Month')
        ax.set_ylabel('Virus Type')
        fig.savefig(os.path.join(static_graph_dir, "aylara_gore_virus_dagilimi.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        report_content += """
        <h3>Monthly Virus Distribution (All Years Total)</h3>
        <img src="static_graphs/aylara_gore_virus_dagilimi.png" alt="Monthly Virus Distribution">
        """

        fig, ax = plt.subplots(figsize=(14, 8))
        year_virus = dataframe.groupby('Year')[virus_cols].sum().T
        sns.heatmap(year_virus, annot=True, fmt='d', cmap='plasma', ax=ax)
        ax.set_title('Yearly Virus Distribution')
        ax.set_xlabel('Year')
        ax.set_ylabel('Virus Type')
        fig.savefig(os.path.join(static_graph_dir, "yillara_gore_virus_dagilimi.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        report_content += """
        <h3>Yearly Virus Distribution</h3>
        <img src="static_graphs/yillara_gore_virus_dagilimi.png" alt="Yearly Virus Distribution">
        """

        if 'Week_of_Year' in dataframe.columns:
            fig, ax = plt.subplots(figsize=(16, 8))
            week_virus = dataframe.groupby('Week_of_Year')[virus_cols].sum().T
            sns.heatmap(week_virus, annot=False, fmt='d', cmap='magma', cbar_kws={'label': 'Case Count'}, ax=ax)
            ax.set_title('Weekly Virus Distribution (All Years Total)')
            ax.set_xlabel('Week of Year')
            ax.set_ylabel('Virus Type')
            fig.savefig(os.path.join(static_graph_dir, "haftalara_gore_virus_dagilimi.png"), bbox_inches='tight', dpi=300)
            plt.close(fig)
            report_content += """
            <h3>Weekly Virus Distribution (All Years Total)</h3>
            <img src="static_graphs/haftalara_gore_virus_dagilimi.png" alt="Weekly Virus Distribution">
            """

        monthly_total_virus = dataframe.set_index('Analyzed_Date')['Total_Virus_Count'].resample('M').sum()
        if not monthly_total_virus.empty:
            fig, ax = plt.subplots(figsize=(15, 7))
            monthly_total_virus.plot(kind='line', marker='o', linestyle='-', color='blue', ax=ax)
            ax.set_title('Monthly Total Virus Case Trend')
            ax.set_xlabel('Date')
            ax.set_ylabel('Total Case Count')
            ax.grid(True)
            fig.savefig(os.path.join(static_graph_dir, "aylik_toplam_virus_trend.png"), bbox_inches='tight', dpi=300)
            plt.close(fig)
            report_content += """
            <h3>Monthly Total Virus Case Trend</h3>
            <img src="static_graphs/aylik_toplam_virus_trend.png" alt="Monthly Total Virus Case Trend">
            """


    # Environmental Factors
    if env_temp_col_name in dataframe.columns and env_hum_col_name in dataframe.columns and virus_cols and dataframe[env_temp_col_name].notna().any() and dataframe[env_hum_col_name].notna().any():
        report_content += """
        <h2>4. Environmental Factors Analysis</h2>
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        sns.histplot(dataframe[env_temp_col_name].dropna(), bins=20, kde=True, color='crimson', ax=ax1)
        ax1.set_title('Monthly Temperature Distribution')
        ax1.set_xlabel(f'{env_temp_col_name} (°C)')
        ax1.set_ylabel('Record Count')

        sns.histplot(dataframe[env_hum_col_name].dropna(), bins=20, kde=True, color='royalblue', ax=ax2)
        ax2.set_title('Monthly Humidity Distribution')
        ax2.set_xlabel(f'{env_hum_col_name} (%)')
        ax2.set_ylabel('Record Count')
        plt.tight_layout()
        fig.savefig(os.path.join(static_graph_dir, "sicaklik_nem_dagilimi.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)
        report_content += """
        <h3>Temperature and Humidity Distribution</h3>
        <img src="static_graphs/sicaklik_nem_dagilimi.png" alt="Temperature and Humidity Distribution">
        """

        temp_bins = pd.cut(dataframe[env_temp_col_name], bins=6)
        if not temp_bins.empty:
            temp_virus = dataframe.groupby(temp_bins)[virus_cols].sum().T
            if not temp_virus.empty:
                fig, ax = plt.subplots(figsize=(14, 8))
                sns.heatmap(temp_virus, annot=True, fmt='d', cmap='coolwarm', ax=ax)
                ax.set_title('Virus Distribution by Temperature Ranges')
                ax.set_xlabel('Temperature Range (°C)')
                ax.set_ylabel('Virus Type')
                fig.savefig(os.path.join(static_graph_dir, "sicaklik_virus_dagilimi.png"), bbox_inches='tight', dpi=300)
                plt.close(fig)
                report_content += """
                <h3>Virus Distribution by Temperature Ranges</h3>
                <img src="static_graphs/sicaklik_virus_dagilimi.png" alt="Virus Distribution by Temperature">
                """

        humidity_bins = pd.cut(dataframe[env_hum_col_name], bins=6)
        if not humidity_bins.empty:
            humidity_virus = dataframe.groupby(humidity_bins)[virus_cols].sum().T
            if not humidity_virus.empty:
                fig, ax = plt.subplots(figsize=(14, 8))
                sns.heatmap(humidity_virus, annot=True, fmt='d', cmap='summer', ax=ax)
                ax.set_title('Virus Distribution by Humidity Ranges')
                ax.set_xlabel('Humidity Range (%)')
                ax.set_ylabel('Virus Type')
                fig.savefig(os.path.join(static_graph_dir, "nem_virus_dagilimi.png"), bbox_inches='tight', dpi=300)
                plt.close(fig)
                report_content += """
                <h3>Virus Distribution by Humidity Ranges</h3>
                <img src="static_graphs/nem_virus_dagilimi.png" alt="Virus Distribution by Humidity">
                """

    if virus_cols and 'Total_Virus_Count' in dataframe.columns:
        report_content += """
        <h2>5. Co-infection Analysis</h2>
        """
        co_infected_count = dataframe[dataframe['Total_Virus_Count'] > 1]
        report_content += f"<p>Patients with co-infections: {len(co_infected_count)} ({(len(co_infected_count)/len(dataframe)*100):.2f}%)</p>"

        # Include association rules table if generated
        data_tables_dir = os.path.join(temp_directory, "data_tables")
        os.makedirs(data_tables_dir, exist_ok=True) # Ensure it exists
        rules_path = os.path.join(data_tables_dir, "association_rules.xlsx")
        try:
            from mlxtend.frequent_patterns import apriori, association_rules
            virus_binary = dataframe[virus_cols].applymap(lambda x: 1 if x > 0 else 0)
            virus_binary = virus_binary[virus_binary.sum(axis=1) > 0]
            if not virus_binary.empty:
                # Use a default min_support for report generation if not explicitly passed
                frequent_itemsets = apriori(virus_binary, min_support=0.005, use_colnames=True)
                if not frequent_itemsets.empty:
                    rules = association_rules(frequent_itemsets, metric="lift", min_threshold=1)
                    rules.sort_values('lift', ascending=False).head(10).to_excel(rules_path, index=False)
                    report_content += """
                    <h3>Virus Association Rules (Top 10 by Lift)</h3>
                    <p>See "data_tables/association_rules.xlsx" for full details.</p>
                    """
        except ImportError:
            report_content += "<p><i>mlxtend library not installed. Association rules could not be generated.</i></p>"
        except Exception as e:
            report_content += f"<p><i>An error occurred during association rule mining for report: {e}</i></p>"

    # --- Report section for Advanced Statistical Analysis (modified) ---
    report_content += """
    <h2>6. Advanced Statistical Analysis</h2>
    <p>Interactive statistical analysis was performed for individual viruses against demographic (Gender, Age Group, Age) and environmental factors (Temperature, Humidity), as well as time-based groupings (Month, Season).</p>
    <p>Results for Chi-square tests (for categorical grouping variables) and Independent t-tests (for continuous grouping variables) were generated interactively in the application, along with relevant visualizations and interpretations.</p>
    """

    numerical_cols = dataframe.select_dtypes(include=np.number).columns.tolist()
    all_numerical_for_clustering_report = [col for col in numerical_cols if col not in ['Total_Virus_Count', 'Year', 'Month', 'Week_of_Year']]

    if all_numerical_for_clustering_report and len(dataframe) > 0:
        report_content += """
        <h2>7. Clustering Analysis</h2>
        <p>Patients were grouped using clustering algorithms to identify distinct patterns based on their numerical features, including virus presence, age, temperature, and humidity.</p>
        """
        # For report, we will generate K-Means plots as they are generally more straightforward for static representation
        # and assume K-Means was run in the app for the report generation button.
        cluster_data_report_alg = dataframe[all_numerical_for_clustering_report].copy()
        for col in all_numerical_for_clustering_report:
            cluster_data_report_alg[col] = pd.to_numeric(cluster_data_report_alg[col], errors='coerce')
        cluster_data_report_alg.dropna(inplace=True)

        if not cluster_data_report_alg.empty and len(cluster_data_report_alg) >= 2: # Ensure at least 2 samples for clustering
            # Default to 4 clusters for report if not specified in advanced options
            num_clusters_report = 4
            kmeans_report = KMeans(n_clusters=num_clusters_report, random_state=42, n_init='auto')
            df_for_cluster_plot_report = dataframe.copy()
            df_for_cluster_plot_report.loc[cluster_data_report_alg.index, 'Cluster'] = kmeans_report.fit_predict(cluster_data_report_alg)
            df_for_cluster_plot_report['Cluster'] = df_for_cluster_plot_report['Cluster'].astype(str)

            # Define features for the scatter plot in the report
            report_x_feature = 'Yaş' if 'Yaş' in df_for_cluster_plot_report.columns else (all_numerical_for_clustering_report[0] if len(all_numerical_for_clustering_report) > 0 else None)
            report_y_feature = env_temp_col_name if env_temp_col_name in df_for_cluster_plot_report.columns else (all_numerical_for_clustering_report[1] if len(all_numerical_for_clustering_report) > 1 else None)

            if report_x_feature and report_y_feature and report_x_feature in df_for_cluster_plot_report.columns and report_y_feature in df_for_cluster_plot_report.columns:
                fig_scatter = px.scatter(
                    df_for_cluster_plot_report.dropna(subset=['Cluster', report_x_feature, report_y_feature]),
                    x=report_x_feature, y=report_y_feature, color='Cluster',
                    title=f'Patient Clusters ({report_x_feature} vs {report_y_feature}) - K-Means',
                    labels={report_x_feature: report_x_feature, report_y_feature: report_y_feature}
                )
                fig_scatter.write_image(os.path.join(static_graph_dir, "cluster_scatter_report.png"))
                report_content += f"""
                <h3>Cluster Visualization ({report_x_feature} vs {report_y_feature})</h3>
                <img src="static_graphs/cluster_scatter_report.png" alt="Cluster Scatter Plot">
                """

            if 'Cluster' in df_for_cluster_plot_report.columns and virus_cols:
                virus_profile_by_cluster_report = df_for_cluster_plot_report.groupby('Cluster')[virus_cols].mean().T
                if not virus_profile_by_cluster_report.empty:
                    fig_virus_heatmap, ax_virus_heatmap = plt.subplots(figsize=(12, 8))
                    sns.heatmap(virus_profile_by_cluster_report, annot=True, fmt=".2f", cmap='Greens', ax=ax_virus_heatmap)
                    ax_virus_heatmap.set_title('Average Virus Presence in Each Cluster')
                    ax_virus_heatmap.set_xlabel('Cluster')
                    ax_virus_heatmap.set_ylabel('Virus Type')
                    fig_virus_heatmap.savefig(os.path.join(static_graph_dir, "cluster_virus_heatmap.png"), bbox_inches='tight', dpi=300)
                    plt.close(fig_virus_heatmap)
                    report_content += """
                    <h3>Average Virus Presence by Cluster</h3>
                    <img src="static_graphs/cluster_virus_heatmap.png" alt="Virus Profile Heatmap">
                    """

            characteristics_features_report = ['Yaş', 'Cinsiyet', env_temp_col_name, env_hum_col_name, 'Month', 'Season']
            available_characteristics_report = [f for f in characteristics_features_report if f in df_for_cluster_plot_report.columns and df_for_cluster_plot_report[f].notna().any()]

            if available_characteristics_report:
                report_content += """<h3>Cluster Characteristics by Other Features</h3>"""
                for feature in available_characteristics_report:
                    if feature in ['Yaş', env_temp_col_name, env_hum_col_name]: # Continuous
                        fig_box, ax_box = plt.subplots(figsize=(10, 6))
                        sns.boxplot(data=df_for_cluster_plot_report, x='Cluster', y=feature, palette='pastel', ax=ax_box)
                        ax_box.set_title(f'{feature} Distribution by Cluster')
                        box_plot_path = os.path.join(static_graph_dir, f"cluster_{feature.replace(' ', '_').replace('/', '_')}_boxplot.png")
                        fig_box.savefig(box_plot_path, bbox_inches='tight', dpi=300)
                        plt.close(fig_box)
                        report_content += f"""
                        <h4>{feature} Distribution</h4>
                        <img src="static_graphs/cluster_{feature.replace(' ', '_').replace('/', '_')}_boxplot.png" alt="{feature} Box Plot">
                        """
                    elif feature in ['Cinsiyet', 'Month', 'Season']: # Categorical
                        fig_count, ax_count = plt.subplots(figsize=(10, 6))
                        sns.countplot(data=df_for_cluster_plot_report, x='Cluster', hue=feature, palette='viridis', ax=ax_count)
                        ax_count.set_title(f'{feature} Distribution by Cluster')
                        ax_count.tick_params(axis='x', rotation=45)
                        count_plot_path = os.path.join(static_graph_dir, f"cluster_{feature.replace(' ', '_').replace('/', '_')}_countplot.png")
                        fig_count.savefig(count_plot_path, bbox_inches='tight', dpi=300)
                        plt.close(fig_count)
                        report_content += f"""
                        <h4>{feature} Distribution</h4>
                        <img src="static_graphs/cluster_{feature.replace(' ', '_').replace('/', '_')}_countplot.png" alt="{feature} Count Plot">
                        """
            else:
                report_content += "<p>No other characteristics available for detailed cluster distribution plots in the report.</p>"
        else:
            report_content += "<p>Insufficient data for clustering analysis in the report.</p>"


    report_content += """
    </body>
    </html>
    """

    report_file_path = os.path.join(temp_directory, "report.html")
    with open(report_file_path, "w", encoding="utf-8") as f:
        f.write(report_content)
    return report_file_path


# Download report
if st.button("Generate Full Report"):
    with st.spinner("Generating report..."):
        # Create a temporary directory for report assets
        with tempfile.TemporaryDirectory() as temp_dir:
            # Generate the HTML report and save all associated plots
            report_path = generate_html_report_streamlit(df, temp_dir, actual_virus_columns, labels, env_temp_col, env_hum_col)

            # Create a zip file
            import zipfile
            zip_path = os.path.join(temp_dir, "report.zip")
            with zipfile.ZipFile(zip_path, 'w') as zipf:
                for root, dirs, files in os.walk(temp_dir):
                    for file in files:
                        # Ensure paths are relative to the zip file's root
                        arcname = os.path.relpath(os.path.join(root, file), temp_dir)
                        zipf.write(os.path.join(root, file), arcname)

            # Provide download link
            with open(zip_path, 'rb') as f:
                bytes_data = f.read()
                b64 = base64.b64encode(bytes_data).decode()
                href = f'<a href="data:application/zip;base64,{b64}" download="respiratory_virus_report.zip">Download Full Report (ZIP)</a>'
                st.markdown(href, unsafe_allow_html=True)
            st.success("Report generated successfully!")

