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
from datetime import datetime, timedelta, date # datetime ve timedelta import edildi
import os
import base64 # Import base64 for the download button
import tempfile
import shutil # Import shutil for directory operations
import statsmodels.api as sm

import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps # Pearsonr için
from sklearn.metrics import roc_curve, roc_auc_score
from statsmodels.stats.multitest import multipletests # Bonferroni düzeltmesi için
from rolling_window_tmp_gem import perform_lagged_correlation_analysis_and_plot

import calendar # Haftalık tarihleri doğru bulmak için gerekli
import io # Bu satırı dosyanızın en başına ekleyin

# tqdm kaldırıldı, Streamlit'in kendi progress barını kullanacağız.

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

def iso_week_start(yil, hafta):
    try:
        return date.fromisocalendar(int(yil), int(hafta), 1)
    except:
        return pd.NaT

def format_p(p):
    try:
        p = float(p)
        if p < 1e-4 and p > 0:
            return "<0.0001"
        elif p == 0:
            return "<1e-16"
        else:
            return f"{p:.4g}"
    except:
        return str(p)

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
            'Adenovirus', 'Coronavirus_HKU1', 'Human_Bocavirus', 'Human_Coronavirus_229E', 'Human_Coronavirus_NL63',
            'Human_Coronavirus_OC43', 'Human_Metapneumovirus',
            'Human_parechovirus', 'Influenza_A', 'Influenza_B',
            'Parainfluenza_Virus_1', 'Parainfluenza_Virus_2',
            'Parainfluenza_Virus_3', 'Parainfluenza_Virus_4',
            'Respiratuvar_sinsityal_virüs_A_B', 'Enterovirus_Rhinovirus', 'SARS-COV-2'
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


# ... (Streamlit uygulamanızın başlangıç kısmı, set_page_config, CSS,
#       load_data fonksiyonu, dosya yükleyici ve df DataFrame'inin oluşturulduğu yer) ...

# df DataFrame'inin yüklendiğinden ve boş olmadığından emin olun.
# Bu kod bloğu, Streamlit'in sidebar'ında dosya yüklendikten sonra çalışmalıdır.
# ... (Streamlit uygulamanızın başlangıç kısmı, set_page_config, CSS,
#       load_data fonksiyonu, dosya yükleyici ve df DataFrame'inin oluşturulduğu yer) ...

# ... (Streamlit uygulamanızın başlangıç kısmı, set_page_config, CSS,
#       load_data fonksiyonu, dosya yükleyici ve df DataFrame'inin oluşturulduğu yer) ...

# df DataFrame'inin yüklendiğinden ve boş olmadığından emin olun.
# Bu kod bloğu, Streamlit'in sidebar'ında dosya yüklendikten sonra çalışmalıdır.
# ... (Streamlit uygulamanızın başlangıç kısmı, set_page_config, CSS,
#       load_data fonksiyonu, dosya yükleyici ve df DataFrame'inin oluşturulduğu yer) ...

# df DataFrame'inin yüklendiğinden ve boş olmadığından emin olun.
# Bu kod bloğu, Streamlit'in sidebar'ında dosya yüklendikten sonra çalışmalıdır.
if 'df' not in locals() or df.empty:
    st.info("Lojistik regresyon analizi yapabilmek için lütfen sol kenar çubuğundan bir Excel dosyası yükleyin.")
else:
    # --- 7. HEADER: Lojistik Regresyon Analizi Bölümü ---
    st.markdown("---")
    st.header("📈 Lojistik Regresyon Analizi")
    st.write("Bu bölüm, çevresel faktörler (nem ve sıcaklık) ile virüs varlığı arasındaki ilişkileri, seçtiğiniz filtrelemelere göre lojistik regresyon kullanarak inceler.")

    st.subheader("Veri Filtreleme Seçenekleri")

    # Veriyi filtrelemek için filtre seçenekleri
    filtered_df = df.copy()

    # Cinsiyet Filtresi (Sütun adı 'Cinsiyet', değerler 'Erkek' ve 'Bayan')
    if 'Cinsiyet' in filtered_df.columns:
        unique_genders = ['Tümü'] + filtered_df['Cinsiyet'].dropna().unique().tolist()
        selected_gender = st.radio(
            "Cinsiyete Göre Filtrele:",
            options=unique_genders,
            key="lr_gender_filter_main"
        )
        
        if selected_gender != 'Tümü':
            filtered_df = filtered_df[filtered_df['Cinsiyet'] == selected_gender]
            st.info(f"Veri seti sadece '{selected_gender}' cinsiyetine sahip hastaları içerecek şekilde filtrelendi. Kalan örnek sayısı: **{len(filtered_df)}**")
    else:
        st.info("Veri setinde 'Cinsiyet' kolonu bulunamadı. Cinsiyet filtresi uygulanamıyor.")


    # Yaş Grubu Filtresi (Sütun adı 'Yaş')
    if 'Yaş' in filtered_df.columns and pd.api.types.is_numeric_dtype(filtered_df['Yaş']):
        min_age_val = int(filtered_df['Yaş'].min()) if not filtered_df['Yaş'].empty else 0
        max_age_val = int(filtered_df['Yaş'].max()) if not filtered_df['Yaş'].empty else 100
        
        age_range = st.slider(
            "Yaş Aralığına Göre Filtrele:",
            min_value=min_age_val,
            max_value=max_age_val,
            value=(min_age_val, max_age_val),
            key="lr_age_filter_main"
        )
        filtered_df = filtered_df[(filtered_df['Yaş'] >= age_range[0]) & (filtered_df['Yaş'] <= age_range[1])]
        st.info(f"Veri seti {age_range[0]}-{age_range[1]} yaş aralığına göre filtrelendi. Kalan örnek sayısı: **{len(filtered_df)}**")
    else:
        st.info("Veri setinde 'Yaş' kolonu bulunamadı veya sayısal değil. Yaş filtresi uygulanamıyor.")
    
    st.write(f"Filtreleme sonrası analiz için kalan toplam örnek sayısı: **{len(filtered_df)}**")

    if filtered_df.empty:
        st.warning("Uygulanan filtreler sonucunda veri seti boş kaldı. Lütfen filtreleme seçimlerinizi gözden geçirin.")
        st.stop()

    st.subheader("Lojistik Regresyon Parametre Seçimi")

    numeric_cols = filtered_df.select_dtypes(include=np.number).columns.tolist()

    if not numeric_cols:
        st.warning("Filtrelenmiş veri setinde sayısal kolon bulunamadı. Nem ve Sıcaklık kolonları seçilemiyor.")
        st.stop()

    env_col1_options = [col for col in numeric_cols if 'nem' in col.lower() or 'humidity' in col.lower()]
    default_nem = env_col1_options[0] if env_col1_options else (numeric_cols[0] if numeric_cols else None)

    env_col2_options = [col for col in numeric_cols if 'sıcaklık' in col.lower() or 'temp' in col.lower()]
    default_sicaklik = env_col2_options[0] if env_col2_options else None
    
    if default_nem and default_nem in env_col2_options:
        temp_env_col2_options = [col for col in env_col2_options if col != default_nem]
        env_col2_options = temp_env_col2_options
        if default_sicaklik == default_nem and env_col2_options:
            default_sicaklik = env_col2_options[0]

    if default_sicaklik is None and len(numeric_cols) > 0:
        if len(numeric_cols) > 1 and nem_kolon == numeric_cols[0]:
            default_sicaklik = numeric_cols[1]
        elif len(numeric_cols) > 0:
            default_sicaklik = numeric_cols[0]
        else:
            default_sicaklik = None


    nem_kolon = st.selectbox(
        "Nem Kolonunu Seçin:",
        options=numeric_cols,
        index=numeric_cols.index(default_nem) if default_nem in numeric_cols else 0,
        key="nem_lr_select"
    )
    sicaklik_kolon = st.selectbox(
        "Sıcaklık Kolonunu Seçin:",
        options=numeric_cols,
        index=numeric_cols.index(default_sicaklik) if default_sicaklik in numeric_cols else (1 if len(numeric_cols)>1 else 0),
        key="sicaklik_lr_select"
    )
    
    if st.button("Lojistik Regresyon Analizi Başlat", key="start_lr_analysis"):
        with st.spinner("Analizler yapılıyor... Bu biraz zaman alabilir."):
            excluded_from_virus_detection = [nem_kolon, sicaklik_kolon, 'Yaş', 'Cinsiyet', 'Month', 'Hafta', 'Yıl', 'Örnek',
                                             'Haftalık_Sıcaklık', 'Dogum Tar.', 'Çalışma Ayı', 'Çalışma Haftası', 'Çalışma Yılı',
                                             ]
            excluded_from_virus_detection = [col for col in excluded_from_virus_detection if col in filtered_df.columns]

            tum_virusler = [col for col in filtered_df.columns
                            if ((filtered_df[col].dropna().isin([0, 1]).all() and not filtered_df[col].dropna().empty) or (filtered_df[col].dropna().unique().tolist() == [0, 1]) or (filtered_df[col].dropna().unique().tolist() == [1, 0]))
                            and col not in excluded_from_virus_detection]

            st.write(f"Tespit edilen tüm potansiyel virüs kolonları (0/1 ikili değer içeren): **{len(tum_virusler)}**")
            if tum_virusler:
                with st.expander("Tüm Virüs Kolonlarını Göster"):
                    st.write(tum_virusler)
            else:
                st.warning("0/1 ikili değer içeren virüs kolonu tespit edilemedi. Lütfen veri setinizi kontrol edin veya sütun adlarını ayarlayın.")
                st.stop()


            virus_cols = [col for col in tum_virusler
                          if filtered_df[col].nunique(dropna=True) == 2]


            st.write(f"Analize alınacak virüs kolonları (pozitif ve negatif gözlemler içeren): **{len(virus_cols)}**")
            if virus_cols:
                with st.expander("Analize Alınacak Virüs Kolonlarını Göster"):
                    st.write(virus_cols)
            else:
                st.warning("Filtreleme sonrası analize uygun (hem 0 hem de 1 içeren) virüs kolonu bulunamadı.")
                st.stop()

            results = []

            progress_bar = st.progress(0)
            status_text = st.empty()

            for i, virus in enumerate(virus_cols):
                status_text.text(f"Analiz ediliyor: **{virus}** ({i+1}/{len(virus_cols)})")
                
                current_model_df = filtered_df[[virus, nem_kolon, sicaklik_kolon]].copy().dropna()
                
                if current_model_df.empty:
                    st.warning(f"Uyarı: '{virus}' için ilgili kolonlarda yeterli (NaN olmayan) veri bulunamadı, bu virüs atlanıyor.")
                    continue
                
                if current_model_df[virus].nunique() < 2:
                    st.warning(f"Uyarı: '{virus}' sütunu sadece tek bir değer içerdiği için (hepsi 0 veya hepsi 1) lojistik regresyon uygulanamıyor. Atlanıyor.")
                    continue

                q_sicaklik = f"Q('{sicaklik_kolon}')"
                q_nem = f"Q('{nem_kolon}')"

                # 1) Sadece sıcaklık modeli
                try:
                    formula = f"{virus} ~ {q_sicaklik}"
                    model = smf.logit(formula=formula, data=current_model_df).fit(disp=0)
                    conf_int = np.exp(model.conf_int()) # Katsayıların CI'larını exponentiate et
                    
                    for param in model.params.index:
                        if param == 'Intercept': continue
                        or_value = np.exp(model.params[param])
                        pval = model.pvalues[param]
                        # Güven aralığı alt ve üst sınırlarını ekliyoruz
                        or_lower_ci = conf_int.loc[param, 0]
                        or_upper_ci = conf_int.loc[param, 1]
                        
                        results.append({
                            'Target_Virus': virus,
                            'Model': 'Sadece Sıcaklık',
                            'Parameter': param,
                            'OR': or_value,
                            'OR_lower_CI': or_lower_ci, # Yeni eklendi
                            'OR_upper_CI': or_upper_ci, # Yeni eklendi
                            'p-value': pval,
                            'N_Samples': len(current_model_df)
                        })
                except Exception as e:
                    pass

                # 2) Sadece nem modeli
                try:
                    formula = f"{virus} ~ {q_nem}"
                    model = smf.logit(formula=formula, data=current_model_df).fit(disp=0)
                    conf_int = np.exp(model.conf_int()) # Katsayıların CI'larını exponentiate et

                    for param in model.params.index:
                        if param == 'Intercept': continue
                        or_value = np.exp(model.params[param])
                        pval = model.pvalues[param]
                        # Güven aralığı alt ve üst sınırlarını ekliyoruz
                        or_lower_ci = conf_int.loc[param, 0]
                        or_upper_ci = conf_int.loc[param, 1]

                        results.append({
                            'Target_Virus': virus,
                            'Model': 'Sadece Nem',
                            'Parameter': param,
                            'OR': or_value,
                            'OR_lower_CI': or_lower_ci, # Yeni eklendi
                            'OR_upper_CI': or_upper_ci, # Yeni eklendi
                            'p-value': pval,
                            'N_Samples': len(current_model_df)
                        })
                except Exception as e:
                    pass

                # 3) Kombinasyon (Sıcaklık ve Nem) modeli
                try:
                    formula = f"{virus} ~ {q_sicaklik} + {q_nem}"
                    model = smf.logit(formula=formula, data=current_model_df).fit(disp=0)
                    conf_int = np.exp(model.conf_int()) # Katsayıların CI'larını exponentiate et

                    for param in model.params.index:
                        if param == 'Intercept': continue
                        or_value = np.exp(model.params[param])
                        pval = model.pvalues[param]
                        # Güven aralığı alt ve üst sınırlarını ekliyoruz
                        or_lower_ci = conf_int.loc[param, 0]
                        or_upper_ci = conf_int.loc[param, 1]

                        results.append({
                            'Target_Virus': virus,
                            'Model': 'Sıcaklık + Nem',
                            'Parameter': param,
                            'OR': or_value,
                            'OR_lower_CI': or_lower_ci, # Yeni eklendi
                            'OR_upper_CI': or_upper_ci, # Yeni eklendi
                            'p-value': pval,
                            'N_Samples': len(current_model_df)
                        })
                except Exception as e:
                    pass
                
                progress_bar.progress((i + 1) / len(virus_cols))
            
            status_text.text("Analiz tamamlandı!")
            progress_bar.empty()

            results_df = pd.DataFrame(results)

            st.subheader("Lojistik Regresyon Sonuçları")
            if not results_df.empty:
                results_df['p-value'] = pd.to_numeric(results_df['p-value'], errors='coerce')
                results_df.dropna(subset=['p-value'], inplace=True) 

                if not results_df.empty:
                    # adj_p-value hesaplaması
                    reject, pvals_corrected, _, _ = multipletests(results_df['p-value'], method='bonferroni', alpha=0.05)
                    results_df['adj_p-value'] = pvals_corrected
                    
                    # Sadece anlamlı olanlar ve artık CI sütunları da dahil
                    significant_df = results_df[results_df['adj_p-value'] < 0.05].copy()
                    
                    # Görüntülenecek sütunları düzenleyebiliriz, örneğin OR_lower_CI ve OR_upper_CI'ı OR'ın yanına getirebiliriz
                    if not significant_df.empty:
                        # Sütun sırasını belirleyelim
                        display_cols = [
                            'Target_Virus', 'Model', 'Parameter', 'OR', 
                            'OR_lower_CI', 'OR_upper_CI', 'p-value', 'adj_p-value', 'N_Samples'
                        ]
                        
                        st.write(f"Bonferroni düzeltmesi sonrası toplam anlamlı regresyon sonucu: **{len(significant_df)}** satır")
                        st.dataframe(significant_df[display_cols].sort_values(by='adj_p-value'))
                        
                        import io
                        output = io.BytesIO()
                        with pd.ExcelWriter(output, engine='openpyxl') as writer:
                            significant_df[display_cols].to_excel(writer, index=False, sheet_name='Anlamlı Sonuçlar')
                        xlsx_data = output.getvalue()


                        excel_buffer = io.BytesIO()
                        significant_yearly_corr_df.to_excel(excel_buffer, index=False)
                        excel_buffer.seek(0) # Buffer'ı başa sarın ki içeriği okunabilsin

                        st.download_button(
                            label="Önemli Yıllık Korelasyonu İndir",
                            data=excel_buffer.getvalue(), # Buffer'ın içeriğini (baytları) data parametresine verin
                            file_name="significant_yearly_correlation.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )
                    else:
                        st.info("Analiz yapıldı fakat Bonferroni düzeltmesi sonrası 0.05'ten küçük p-değerine sahip anlamlı sonuç bulunamadı.")
                else:
                    st.warning("p-değerleri dönüştürüldükten sonra sonuç DataFrame'i boş kaldı. Lütfen verilerinizi kontrol edin.")
            else:
                st.warning("Hiçbir lojistik regresyon modeli oluşturulamadı. Lütfen giriş verilerinizi, seçilen parametreleri ve filtreleri kontrol edin.")

# ... (Bu noktanın altında genellikle generate_html_report_streamlit fonksiyon tanımı ve
#       if st.button("Generate Full Report"): bloğu yer alır.) ...


# --- YENİ BÖLÜM: Gecikmeli ve Dönemsel Korelasyon Analizi ---
st.markdown("---")
st.header("⏳ Gecikmeli ve Dönemsel Korelasyon Analizi")
st.write("Bu bölüm, çevresel faktörler ile virüs prevalansı arasındaki zaman gecikmeli ilişkileri ve bu ilişkilerin zaman içinde nasıl değiştiğini inceler.")


# Make a copy to avoid modifying the main df directly for this section's pre-processing
analysis_df = df.copy()

# Ensure 'Yıl' and 'Hafta' are available and correct type for aggregation
if 'Yıl' not in analysis_df.columns:
    analysis_df['Yıl'] = analysis_df['Çalisma Tar.'].dt.year
if 'Hafta' not in analysis_df.columns:
    analysis_df['Hafta'] = analysis_df['Çalisma Tar.'].dt.isocalendar().week.astype(int)

# Select Environmental Column
env_cols = analysis_df.select_dtypes(include=np.number).columns.tolist()
# Filter out date-related numeric cols if they aren't true environmental factors
env_cols = [col for col in env_cols if col not in ['Yıl', 'Hafta', 'Ay', 'Örnek', 'Dogum Tar.']] # 'Dogum Tar.' added if it's numeric

# Try to suggest a default environmental column
default_env_col = None
if 'Haftalık_Sıcaklık' in env_cols:
    default_env_col = 'Haftalık_Sıcaklık'
elif 'Temp_ay_ortalaması' in env_cols:
    default_env_col = 'Temp_ay_ortalaması'
elif 'Sicaklik' in env_cols:
    default_env_col = 'Sicaklik'
elif env_cols:
    default_env_col = env_cols[0]

if not default_env_col:
    st.warning("Analiz için uygun sayısal çevresel faktör kolonu bulunamadı. Lütfen veri setinizi kontrol edin.")
    # Removed st.stop() to allow the app to continue loading
    # st.stop()
else:
    selected_env_col = st.selectbox(
        "Çevresel Faktör Kolonunu Seçin:",
        options=env_cols,
        index=env_cols.index(default_env_col) if default_env_col in env_cols else 0,
        key="lag_env_select"
    )
    
    st.subheader("Parametre Seçimi")
    
    # Determine virus columns for this section (re-use actual_virus_columns from load_data)
    # lag_virus_cols = actual_virus_columns # This variable needs to be defined in the context where this snippet is used. Assuming it's defined elsewhere.
    lag_virus_cols = actual_virus_columns # Placeholder, replace with actual_virus_columns from your full app

    if not lag_virus_cols:
        st.warning("0/1 ikili değer içeren virüs kolonu tespit edilemedi. Lütfen veri setinizi kontrol edin veya sütun adlarını ayarlayın.")
        # st.stop()
        
    # UI for Lagged Correlation
    max_lag_weeks = st.slider("Maksimum Gecikme Haftası (0-8):", min_value=0, max_value=12, value=8, key="max_lag_weeks")

    # UI for Rolling Window
    st.markdown("---")
    st.subheader("Dönemsel (Rolling Window) Korelasyon Ayarları")
    


    # User selects the specific lag to use for rolling window analysis
    available_lags = list(range(max_lag_weeks + 1))
    
    if st.button("Korelasyon Analizlerini Başlat", key="start_correlation_analysis"):
        with st.spinner("Korelasyon analizleri yapılıyor ve grafikler oluşturuluyor..."):
            # --- Haftalık Veri Agregasyonu ---
            # Bu kısım, yıllık ve dönemsel analizler için temel veriyi toplar
            weekly_data_overall = analysis_df.groupby(['Yıl', 'Hafta']).agg(
                Total_Samples=('Örnek', 'count'), 
                **{f'Positive_{v}': (v, 'sum') for v in lag_virus_cols},
                **{f'Avg_{selected_env_col}': (selected_env_col, 'mean')} 
            ).reset_index()

            weekly_data_overall['Hafta_Başı_Tarihi'] = weekly_data_overall.apply(
                lambda x: iso_week_start(x['Yıl'], x['Hafta']), axis=1
            )

            for v in lag_virus_cols:
                weekly_data_overall[f'Prevalence_{v}'] = weekly_data_overall[f'Positive_{v}'] / weekly_data_overall['Total_Samples']
            
            # Eksik haftalık verileri düzgünleştirmek için
            # Eğer haftalık_veri_overall boşsa veya gerekli sütunlar yoksa uyarı ver
            if weekly_data_overall.empty or f'Avg_{selected_env_col}' not in weekly_data_overall.columns:
                st.warning("Haftalık veri toplama başarısız oldu veya çevresel kolon bulunamadı. Lütfen veri setinizi ve kolon adlarını kontrol edin.")
                # st.stop() # Hatanın devamını engellemek için durdurma
            else:
                # Create a date for each week for plotting and rolling window indexing
                # Hafta sonu tarihlerini hesaplamak için calendar modülünü kullanabiliriz.
                # datetime.fromisocalendar, ISO yıl, hafta ve haftanın günü ile tarih oluşturur.
                # Pazartesiyi (1) kullanarak haftanın başlangıcını bulalım.
                # deneme_guncel_tmp2.py (yaklaşık olarak 1018. satır - bu kısmı komple değiştirin)

                # Yıl ve hafta numarasını birleştirerek pd.to_datetime için uygun bir string formatı oluşturun
                weekly_data_overall['YearWeekDay_Str'] = weekly_data_overall['Yıl'].astype(str) + '-W' + \
                                                        weekly_data_overall['Hafta'].astype(str).str.zfill(2) + '-1'

                # Oluşturulan string'i tarihe çevirin. '%Y-W%W-%w' formatı ISO yıl, hafta ve hafta içi (Pazartesi için 1) anlamına gelir.
                # errors='coerce' sayesinde geçersiz haftalar NaT (Not a Time) olarak işaretlenir.
                weekly_data_overall['Week_Start_Date'] = pd.to_datetime(weekly_data_overall['YearWeekDay_Str'], format='%Y-W%W-%w', errors='coerce')

                # Geçici olarak oluşturduğumuz 'YearWeekDay_Str' sütununu kaldırın
                weekly_data_overall.drop(columns=['YearWeekDay_Str'], inplace=True)

                # Geçersiz veya dönüştürülemeyen tarihler içeren satırları kaldırın (NaT değerleri)
                weekly_data_overall.dropna(subset=['Week_Start_Date'], inplace=True)
                    
                    
                
                weekly_data_overall.sort_values(by='Week_Start_Date', inplace=True)
                weekly_data_overall['Week_End_Date'] = weekly_data_overall['Week_Start_Date'] + timedelta(days=6)


                st.subheader("Yıllık Gecikmeli Korelasyon Sonuçları")
                # --- Yıllık Gecikmeli Korelasyon Analizi ---
                yearly_correlation_results = []
                years = weekly_data_overall['Yıl'].unique()

                for year in years:
                    yearly_df = weekly_data_overall[weekly_data_overall['Yıl'] == year].copy()
                    
                    if yearly_df.empty: continue

                    for virus_col in lag_virus_cols:
                        # Virüs prevalans kolonu yoksa atla
                        if f'Prevalence_{virus_col}' not in yearly_df.columns:
                            continue

                        for lag in range(max_lag_weeks + 1):
                            if len(yearly_df) > lag:
                                lagged_env_col_data = yearly_df[f'Avg_{selected_env_col}'].shift(lag)
                                
                                temp_df = pd.DataFrame({
                                    'Virüs_Haftası': yearly_df['Hafta'],
                                    'Virüs_Tarihi': yearly_df['Hafta_Başı_Tarihi'],
                                    'Çevresel_Haftası': yearly_df['Hafta'] - lag,
                                    'Çevresel_Tarihi': yearly_df['Hafta_Başı_Tarihi'].shift(lag),
                                    'Prevalence': yearly_df[f'Prevalence_{virus_col}'],
                                    'Lagged_Env': lagged_env_col_data
                                }).dropna()

                                # Korelasyon için en az 2 veri noktası ve her iki seride de varyasyon olmalı
                                if len(temp_df) > 1 and temp_df['Prevalence'].nunique() > 1 and temp_df['Lagged_Env'].nunique() > 1:
                                    try:
                                        r_value, p_value = sps.pearsonr(temp_df['Lagged_Env'], temp_df['Prevalence'])
                                        
                                        sig_level = ''
                                        if p_value < 0.001: sig_level = '***'
                                        elif p_value < 0.01: sig_level = '**'
                                        elif p_value < 0.05: sig_level = '*'

                                        yearly_correlation_results.append({
                                            'Virüs': virus_col,
                                            'Yıl': year,
                                            'Gecikme (Hafta)': lag,
                                            f'{selected_env_col}_Pearson_R': r_value,
                                            'p-değeri': p_value,
                                            'Anlamlılık': sig_level,
                                            'Örnek_Sayısı': len(temp_df),
                                            'Başlangıç_Tarihi': temp_df['Virüs_Tarihi'].min(),
                                            'Bitiş_Tarihi': temp_df['Virüs_Tarihi'].max(),
                                            'Çevresel_Başlangıç_Tarihi': temp_df['Çevresel_Tarihi'].min(),
                                            'Çevresel_Bitiş_Tarihi': temp_df['Çevresel_Tarihi'].max(),
                                        })
                                    except Exception as e:
                                        pass

                yearly_corr_df = pd.DataFrame(yearly_correlation_results)

                if not yearly_corr_df.empty:
                    if 'p-değeri' in yearly_corr_df.columns and not yearly_corr_df['p-değeri'].empty:
                        # Bonferroni düzeltmesi (Tüm testler için)
                        reject, pvals_corrected, _, _ = multipletests(yearly_corr_df['p-değeri'], method='bonferroni', alpha=0.05)
                        yearly_corr_df['Ayarlı_p-değeri'] = pvals_corrected
                        
                        significant_yearly_corr_df = yearly_corr_df[yearly_corr_df['Ayarlı_p-değeri'] < 0.05].copy()
                        
                        if not significant_yearly_corr_df.empty:
                            st.dataframe(significant_yearly_corr_df.sort_values(by=['Yıl', 'Ayarlı_p-değeri']))
                            
                            # Yıllık Korelasyon sonuçlarını bellekteki bir Excel dosyasına kaydedin
                            yearly_corr_excel_buffer = io.BytesIO()
                            significant_yearly_corr_df.to_excel(yearly_corr_excel_buffer, index=False)
                            yearly_corr_excel_buffer.seek(0) # Buffer'ı başa sarın

                            st.download_button(
                                label="Yıllık Korelasyon Sonuçlarını Excel İndir",
                                data=yearly_corr_excel_buffer.getvalue(), # Bellekteki Excel içeriğini verin
                                file_name="yillik_gecikmeli_korelasyon_sonuclari.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                            )
                        else:
                            st.info("Ayarlı p-değeri < 0.05 olan anlamlı yıllık gecikmeli korelasyon bulunamadı.")
                    else:
                        st.info("Yıllık korelasyon analizi için p-değeri hesaplanamadı.")
                else:
                    st.info("Yıllık gecikmeli korelasyon analizi için yeterli veri bulunamadı veya bir hata oluştu.")


                st.subheader("Dönemsel (Rolling Window) Korelasyon Sonuçları")
                # --- Dönemsel (Rolling Window) Korelasyon Analizi ---
                rolling_window_results = []

                import io
                from scipy.stats import pearsonr
                from statsmodels.stats.multitest import multipletests

                st.subheader("Dönemsel (Rolling Window) Korelasyon Sonuçları (Tüm Virüsler, Tüm Lag/Window)")
                st.info("Sadece Bonferroni düzeltmesi ile anlamlı çıkan korelasyonlar gösterilmektedir (Bonferroni p<0.05).")

                results = []
                for virus_col in lag_virus_cols:
                    for lag in range(2, 5):        # Sadece 2-4 hafta lag!
                        for window in range(3, 6): # Sadece 3-5 hafta window!
                            dfw = weekly_data_overall.sort_values('Week_Start_Date').reset_index(drop=True)
                            lagged_env = dfw[f'Avg_{selected_env_col}'].shift(lag)
                            for i in range(len(dfw) - window - lag + 1):
                                v_window = dfw.iloc[i + lag : i + lag + window]
                                e_window = dfw.iloc[i : i + window]
                                if v_window[f'Prevalence_{virus_col}'].nunique() < 2 or e_window[f'Avg_{selected_env_col}'].nunique() < 2:
                                    continue
                                try:
                                    r, p = pearsonr(e_window[f'Avg_{selected_env_col}'], v_window[f'Prevalence_{virus_col}'])
                                except:
                                    continue
                                if p < 0.05:
                                    results.append({
                                        'Virüs': virus_col,
                                        'Lag': lag,
                                        'Pencere': window,
                                        'Çevresel_Başlangıç': e_window['Week_Start_Date'].iloc[0].date(),
                                        'Çevresel_Bitiş': e_window['Week_Start_Date'].iloc[-1].date(),
                                        'Virüs_Başlangıç': v_window['Week_Start_Date'].iloc[0].date(),
                                        'Virüs_Bitiş': v_window['Week_Start_Date'].iloc[-1].date(),
                                        'Korelasyon': f"{r:.3f}",
                                        'p-değeri': p,
                                        'Çevre_Başlangıç_Degeri': round(e_window[f'Avg_{selected_env_col}'].iloc[0], 2),
                                        'Çevre_Bitiş_Degeri': round(e_window[f'Avg_{selected_env_col}'].iloc[-1], 2),
                                        'Virüs_Başlangıç_Vaka': int(v_window[f'Positive_{virus_col}'].iloc[0]),
                                        'Virüs_Bitiş_Vaka':     int(v_window[f'Positive_{virus_col}'].iloc[-1]),
                                    })
                rolling_sig_df = pd.DataFrame(results)

                if not rolling_sig_df.empty:
                    N = len(rolling_sig_df)
                    # Bonferroni ve FDR düzeltmeleri
                    rolling_sig_df['Bonferroni_p'] = (rolling_sig_df['p-değeri'].astype(float) * N).clip(upper=1)
                    _, fdr_p, _, _ = multipletests(rolling_sig_df['p-değeri'].astype(float), alpha=0.05, method='fdr_bh')
                    rolling_sig_df['FDR_p'] = fdr_p
                    # Bilimsel p-değeri gösterimi
                    for col in ['p-değeri', 'Bonferroni_p', 'FDR_p']:
                        rolling_sig_df[col] = rolling_sig_df[col].apply(format_p)
                    # Sadece Bonferroni anlamlıları tut
                    rolling_sig_df = rolling_sig_df[rolling_sig_df['Bonferroni_p'].apply(lambda x: not x.startswith('>') and float(x.replace('<','').replace('e','E').replace('.','0.')) < 0.05 if x[0].isdigit() else True)]
                    rolling_sig_df['Anlamlılık'] = '***'
                    # Sonuçları göster ve indir
                    if not rolling_sig_df.empty:
                        st.dataframe(rolling_sig_df.sort_values('p-değeri'))
                        import io
                        buffer = io.BytesIO()
                        rolling_sig_df.to_excel(buffer, index=False)
                        st.download_button("Sonuçları Excel Olarak İndir", buffer.getvalue(), "donemsel_korelasyon_bonferroni.xlsx")
                    else:
                        st.warning("Belirtilen lag ve pencere aralığında, Bonferroni düzeltmesiyle anlamlı korelasyon bulunamadı.")
                else:
                    st.warning("Hiç anlamlı korelasyonlu dönem bulunamadı.")




            
    else:
        st.info("Lütfen 'Korelasyon Analizlerini Başlat' butonuna basarak analizleri başlatın.")

    # Sonuçları Yorumlama Rehberi
    st.markdown("---")
    st.subheader("Sonuçları Yorumlama Rehberi")
    st.info("""
        **Pearson R Değeri:**
        * **Pozitif (R > 0):** Çevresel faktör arttıkça virüs prevalansı da artma eğilimindedir. (Örn: Sıcaklık arttıkça virüs artıyor)
        * **Negatif (R < 0):** Çevresel faktör arttıkça virüs prevalansı azalma eğilimindedir. (Örn: Sıcaklık arttıkça virüs azalıyor)
        * **Güç (Mutlak Değer):**
            * 0.0 - 0.2: Çok zayıf/İhmal edilebilir
            * 0.2 - 0.4: Zayıf
            * 0.4 - 0.6: Orta
            * 0.6 - 0.8: Güçlü
            * 0.8 - 1.0: Çok Güçlü
        
        **Gecikme (Hafta):**
        * Gecikme, çevresel faktördeki bir değişimin virüs prevalansını etkilemesinin ne kadar zaman aldığını gösterir. Örneğin, 4 hafta gecikmeli bir korelasyon, 4 hafta önceki sıcaklığın şimdiki virüs prevalansıyla ilişkili olduğunu gösterir.

        **Dönemsel (Rolling Window) Korelasyon:**
        * Bu analiz, belirli bir çevresel faktör-virüs ilişkisinin gücünün ve yönünün zaman içinde (mevsimden mevsime) nasıl değiştiğini gösterir. Aynı ilişkinin yılın bir bölümünde çok güçlü ve anlamlıyken, başka bir bölümünde zayıf veya anlamsız olabileceğini gözlemleyebilirsiniz.
        * Timeline grafiği, anlamlı ilişkilerin hangi zaman aralıklarında ortaya çıktığını görselleştirir.

        **Unutmayın:** Korelasyon nedensellik değildir. Bu analizler, çevresel faktörler ile virüs yayılımı arasındaki istatistiksel ilişkileri ortaya koyar, ancak birinin diğerine doğrudan neden olduğunu kanıtlamaz.
    """)


# ... (deneme_guncel_tmp2.py dosyanızın kalan kodu)
# Satır 400 civarı veya rapor indirme düğmesinden hemen önce:

# Mevcut "Download report" kısmından ÖNCE bu yeni rolling_window analiz bloğunu ekleyin.

# ----------------------------------------------------------------------------------------------------------------------
# | Yeni Eklenen Bölüm: Gecikmeli Korelasyon Analizi (Rolling Window)                                                  |
# ----------------------------------------------------------------------------------------------------------------------

# final_merged_df_weekly yerine st.session_state.df kontrol ediyoruz.
# df'in st.session_state'te saklandığını ve ana işleme sonrası kullanıma hazır olduğunu varsayıyoruz.
if 'df' in st.session_state and st.session_state.df is not None: # 'df' ana DataFrame'inizse bu kontrolü kullanın
    processed_df_for_analysis = st.session_state.df # Artık bu bizim analiz DataFrame'imiz

    st.sidebar.markdown("---")
    st.sidebar.header("Gecikmeli Korelasyon Analizi (Rolling Window)")

    # Virüs sütunlarını dinamik olarak al
    available_virus_columns = [col for col in processed_df_for_analysis.columns if col.endswith('_Prevalansı')]

    if not available_virus_columns:
        st.sidebar.warning("Analiz için uygun virüs prevalansı sütunu bulunamadı (Örn: 'VirüsAdı_Prevalansı').")
        selected_virus_columns = []
    else:
        selected_virus_columns = st.sidebar.multiselect(
            "Analiz Edilecek Virüs(ler)i Seçin:",
            options=available_virus_columns,
            default=available_virus_columns
        )

    rolling_window_size_input = st.sidebar.slider(
        "Dönemsel Pencere Boyutu (Hafta):",
        min_value=0,
        max_value=52,
        value=12,
        step=1
    )
    max_lag_weeks_input = st.sidebar.slider(
        "Maksimum Gecikme (Hafta):",
        min_value=0,
        max_value=20,
        value=8,
        step=1
    )

    perform_rolling_analysis = st.sidebar.button("Gecikmeli Korelasyon Analizini Başlat")

    if perform_rolling_analysis and selected_virus_columns:
        st.subheader("Gecikmeli Korelasyon Analizi Sonuçları")
        st.info("Analiz başlatılıyor, bu işlem biraz zaman alabilir...")

        for virus_col in actual_virus_columns:
            st.subheader(f"Dönemsel Gecikmeli Analiz Sonuçları: {virus_col}")
            overall_corr_df, yearly_corr_df, rolling_corr_df, figures_dict = \
                perform_lagged_correlation_analysis_and_plot(
                    base_output_excel_folder="streamlit_rolling_reports",
                    base_output_plots_folder="streamlit_rolling_plots",
                    max_lag_weeks=max_lag_weeks,
                    rolling_window_size=rolling_window_size,
                    df_input=weekly_data_overall,
                    virus_column_name=virus_col,
                    env_factor_column_name=selected_env_col,
                    save_to_file=False
                )
            if rolling_corr_df is not None and not rolling_corr_df.empty:
                st.dataframe(rolling_corr_df)
                if 'rolling_corr_trend' in figures_dict:
                    st.pyplot(figures_dict['rolling_corr_trend'])
            else:
                st.info(f"{virus_col} için dönemsel analizde yeterli veri bulunamadı.")

            if overall_corr_df is not None:
                st.write(f"**{virus_col} - Genel Gecikmeli Korelasyon Tablosu**")
                st.dataframe(overall_corr_df)
                if 'overall_heatmap' in figures_dict and figures_dict['overall_heatmap'] is not None:
                    st.pyplot(figures_dict['overall_heatmap'])
                    plt.close(figures_dict['overall_heatmap'])

                if yearly_corr_df is not None and not yearly_corr_df.empty:
                    st.write(f"**{virus_col} - Yıllık Gecikmeli Korelasyon Tablosu**")
                    st.dataframe(yearly_corr_df)
                    if 'yearly_heatmap' in figures_dict and figures_dict['yearly_heatmap'] is not None:
                        st.pyplot(figures_dict['yearly_heatmap'])
                        plt.close(figures_dict['yearly_heatmap'])

                    if rolling_corr_df is not None and not rolling_corr_df.empty:
                        st.dataframe(rolling_corr_df.head())
                        # Görsel varsa göster
                        if 'rolling_corr_trend' in figures_dict:
                            st.pyplot(figures_dict['rolling_corr_trend'])


                    st.write(f"**{virus_col} - Haftalık Prevalans ve Sıcaklık Trendleri**")
                    if 'prevalance_temp_trend' in figures_dict and figures_dict['prevalance_temp_trend'] is not None:
                        st.pyplot(figures_dict['prevalance_temp_trend'])
                        plt.close(figures_dict['prevalance_temp_trend'])

                    st.write(f"**{virus_col} - Dönemsel Korelasyon Trendi (Pencere: {rolling_window_size_input} Hafta)**")
                    if 'rolling_corr_trend' in figures_dict and figures_dict['rolling_corr_trend'] is not None:
                        st.pyplot(figures_dict['rolling_corr_trend'])
                        plt.close(figures_dict['rolling_corr_trend'])

                st.write(f"**{virus_col} - Aylık Prevalans Dağılımı**")
                if 'monthly_boxplot' in figures_dict and figures_dict['monthly_boxplot'] is not None:
                    st.pyplot(figures_dict['monthly_boxplot'])
                    plt.close(figures_dict['monthly_boxplot'])

                st.write(f"**{virus_col} - Prevalans vs. Sıcaklık Serpme Grafiği**")
                if 'scatter_plot' in figures_dict and figures_dict['scatter_plot'] is not None:
                    st.pyplot(figures_dict['scatter_plot'])
                    plt.close(figures_dict['scatter_plot'])

            else:
                st.error(f"'{virus_col}' için gecikmeli korelasyon analizi başarısız oldu veya veri bulunamadı.")

        st.success("Tüm gecikmeli korelasyon analizleri tamamlandı!")
    elif perform_rolling_analysis and not selected_virus_columns:
        st.warning("Lütfen analiz edilecek en az bir virüs sütunu seçin.")

else:
    st.info("Gecikmeli Korelasyon Analizi bölümünü görmek için lütfen yukarıdan verilerinizi yükleyin ve işleyin.")

# ----------------------------------------------------------------------------------------------------------------------
# | Yeni Eklenen Bölüm Sonu                                                                                            |
# ----------------------------------------------------------------------------------------------------------------------


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

