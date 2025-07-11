# rolling_window.py

import pandas as pd
import numpy as np
import scipy.stats as sps
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
import shutil
import calendar # calendar eklendi

warnings.filterwarnings('ignore')
pd.set_option('display.max_columns', None)
plt.style.use('ggplot')
sns.set_style("whitegrid")

def sanitize_filename(filename):
    problematic_chars = ['/', '\\', ':', '*', '?', '"', '<', '>', '|']
    for char in problematic_chars:
        filename = filename.replace(char, '_')
    return filename

# perform_lagged_correlation_analysis_and_plot fonksiyonunu aşağıdaki gibi güncelleyin
def perform_lagged_correlation_analysis_and_plot(base_output_excel_folder, base_output_plots_folder,
                                          max_lag_weeks=12, rolling_window_size=12, # max_lag_weeks changed to 12
                                          df_input=None, virus_column_name='Adenovirus_Prevalans',
                                          env_factor_column_name='Haftalık_Sıcaklık',
                                          rolling_analysis_lag=0, # Yeni parametre: rolling_analysis_lag
                                          save_to_file=True): # Yeni parametre: save_to_file
    """
    Birleştirilmiş haftalık veri üzerinde genel ve yıla özel gecikmeli korelasyon analizi yapar,
    sonuçları Excel'e kaydeder ve önemli görselleştirmeleri .png olarak kaydeder.
    'Haftalık_Sıcaklık' sütununu çevresel faktör olarak kullanır.
    Dönemsel (rolling window) korelasyon analizi eklenmiştir.

    Args:
        base_output_excel_folder (str): Korelasyon sonuçlarının ve özet raporların kaydedileceği ana Excel klasörü.
        base_output_plots_folder (str): Görselleştirmelerin kaydedileceği ana klasör yolu.
        max_lag_weeks (int): Analiz edilecek maksimum gecikme süresi (hafta cinsinden).
        rolling_window_size (int): Dönemsel korelasyon analizi için pencere boyutu (hafta cinsinden).
        df_input (pd.DataFrame, optional): Doğrudan analiz edilecek DataFrame. Eğer sağlanırsa merged_excel_path yerine kullanılır.
        virus_column_name (str): Analiz edilecek virüs prevalans sütununun adı.
        save_to_file (bool): Sonuçları dosyalara kaydedip kaydetmeyeceğini belirler. True ise kaydeder.

    Returns:
        tuple: (overall_corr_df, yearly_corr_df, rolling_corr_df, figures_dict)
               figures_dict, matplotlib fig objelerini içerir.
               Hata durumunda (None, None, None, None) döndürür.
    """
    # Excel kaydetme yolunu kaldırdık, sadece df_input kullanacağız
    merged_excel_path = None # Artık kullanılmıyor

    overall_corr_df = None
    yearly_corr_df = None
    rolling_corr_df = None
    figures_dict = {}

    try:
        if df_input is not None:
            print("INFO: Sağlanan DataFrame kullanılıyor.")
            df_merged = df_input.copy() # Orijinal DataFrame'i değiştirmemek için kopya al
        else:
            print("HATA: Analiz için bir DataFrame (df_input) sağlanmalıdır.")
            return (None, None, None, None)
        
        print(f"Birleştirilmiş dosya sütunları: {df_merged.columns.tolist()}")

        # Gerekli sütunların varlığını kontrol et
        required_cols = ['Date', 'Hafta', 'Yıl', 'Haftalık_Sıcaklık', virus_column_name]
        for col in required_cols:
            if col not in df_merged.columns:
                print(f"HATA: Gerekli sütun '{col}' bulunamadı.")
                return (None, None, None, None)

        # Sayısal sütunlardaki virgül (,) ondalık ayırıcıyı nokta (.) ile değiştir
        numerical_cols = ['Haftalık_Sıcaklık', virus_column_name]
        for col in numerical_cols:
            if df_merged[col].dtype == 'object':
                df_merged[col] = df_merged[col].str.replace(',', '.', regex=False).astype(float)
            elif df_merged[col].dtype != 'float64' and df_merged[col].dtype != 'int64':
                df_merged[col] = pd.to_numeric(df_merged[col], errors='coerce')
        
        df_merged.dropna(subset=numerical_cols, inplace=True) # NaN değerleri kaldır
        df_merged.reset_index(drop=True, inplace=True)

        # Çıktı klasörlerini oluştur (sadece save_to_file True ise)
        if save_to_file:
            output_excel_folder = os.path.join(base_output_excel_folder, f"rolling_window_{rolling_window_size}W_{sanitize_filename(virus_column_name)}")
            output_plots_folder = os.path.join(base_output_plots_folder, f"rolling_window_{rolling_window_size}W_{sanitize_filename(virus_column_name)}")

            if os.path.exists(output_excel_folder):
                shutil.rmtree(output_excel_folder)
            os.makedirs(output_excel_folder, exist_ok=True)

            if os.path.exists(output_plots_folder):
                shutil.rmtree(output_plots_folder)
            os.makedirs(output_plots_folder, exist_ok=True)
            print(f"Çıktılar '{output_excel_folder}' ve '{output_plots_folder}' klasörlerine kaydedilecektir.")


        # --- Genel Gecikmeli Korelasyon Analizi ---
        print("\n--- Genel Gecikmeli Korelasyon Analizi Başlatılıyor ---")
        lags = range(max_lag_weeks + 1)
        correlations = []
        p_values = []

        for lag in lags:
            if lag == 0:
                corr, p_val = sps.pearsonr(df_merged[virus_column_name], df_merged['Haftalık_Sıcaklık'])
            else:
                # Gecikmeli seriyi oluştur
                df_lagged = df_merged.copy()
                df_lagged['Lagged_Sıcaklık'] = df_lagged['Haftalık_Sıcaklık'].shift(lag)
                df_lagged.dropna(subset=[virus_column_name, 'Lagged_Sıcaklık'], inplace=True)
                
                if len(df_lagged) > 1: # En az iki veri noktası olmalı
                    corr, p_val = sps.pearsonr(df_lagged[virus_column_name], df_lagged['Lagged_Sıcaklık'])
                else:
                    corr, p_val = np.nan, np.nan # Yeterli veri yoksa NaN
            correlations.append(corr)
            p_values.append(p_val)

        overall_corr_df = pd.DataFrame({
            'Gecikme (Hafta)': lags,
            'Korelasyon (r)': correlations,
            'P-Değeri': p_values
        })
        overall_corr_df['Anlamlı (p<0.05)'] = overall_corr_df['P-Değeri'] < 0.05
        print("Genel Gecikmeli Korelasyon Sonuçları:\n", overall_corr_df)

        if save_to_file:
            overall_corr_df.to_excel(os.path.join(output_excel_folder, f'{sanitize_filename(virus_column_name)}_Overall_Lagged_Correlation.xlsx'), index=False)


        # --- Yıllık Gecikmeli Korelasyon Analizi ---
        print("\n--- Yıllık Gecikmeli Korelasyon Analizi Başlatılıyor ---")
        all_yearly_corrs = []
        for year in sorted(df_merged['Yıl'].unique()):
            df_year = df_merged[df_merged['Yıl'] == year].copy()
            if len(df_year) < max_lag_weeks + 1: # Yeterli veri yoksa atla
                print(f"UYARI: {year} yılı için yeterli veri yok, yıllık analiz atlandı.")
                continue

            yearly_correlations = []
            yearly_p_values = []
            start_weeks = []
            end_weeks = []

            for lag in lags:
                if lag == 0:
                    used_weeks = df_year['Hafta']
                    corr, p_val = sps.pearsonr(df_year[virus_column_name], df_year['Haftalık_Sıcaklık'])
                else:
                    df_lagged_year = df_year.copy()
                    df_lagged_year['Lagged_Sıcaklık'] = df_lagged_year['Haftalık_Sıcaklık'].shift(lag)
                    df_lagged_year.dropna(subset=[virus_column_name, 'Lagged_Sıcaklık'], inplace=True)
                    used_weeks = df_lagged_year['Hafta']
                    if len(df_lagged_year) > 1:
                        corr, p_val = sps.pearsonr(df_lagged_year[virus_column_name], df_lagged_year['Lagged_Sıcaklık'])
                    else:
                        corr, p_val = np.nan, np.nan

                # Yeni eklenenler
                if len(used_weeks) > 0:
                    start_weeks.append(used_weeks.min())
                    end_weeks.append(used_weeks.max())
                else:
                    start_weeks.append(np.nan)
                    end_weeks.append(np.nan)
                yearly_correlations.append(corr)
                yearly_p_values.append(p_val)
            
            yearly_df = pd.DataFrame({
                'Yıl': year,
                'Gecikme (Hafta)': lags,
                'Korelasyon (r)': yearly_correlations,
                'P-Değeri': yearly_p_values,
                'Başlangıç_Haftası': start_weeks,
                'Bitiş_Haftası': end_weeks,
            })
            yearly_df['Anlamlı (p<0.05)'] = yearly_df['P-Değeri'] < 0.05
            all_yearly_corrs.append(yearly_df)

        yearly_corr_df = pd.concat(all_yearly_corrs, ignore_index=True) if all_yearly_corrs else pd.DataFrame()
        if not yearly_corr_df.empty:
            print("Yıllık Gecikmeli Korelasyon Sonuçları (İlk 5 Satır):\n", yearly_corr_df.head())
            if save_to_file:
                yearly_corr_df.to_excel(os.path.join(output_excel_folder, f'{sanitize_filename(virus_column_name)}_Yearly_Lagged_Correlation.xlsx'), index=False)
        else:
            print("UYARI: Yıllık korelasyon analizi için yeterli veri bulunamadı.")

        # --- Dönemsel (Rolling Window) Korelasyon Analizi ---
        print(f"\n--- Dönemsel (Rolling Window) Korelasyon Analizi Başlatılıyor (Pencere Boyutu: {rolling_window_size} Hafta) ---")
        
        # Tarih ve hafta sütunlarını birleştirerek benzersiz haftalık başlangıç noktaları oluştur
        df_merged['YearWeek'] = df_merged['Yıl'].astype(str) + '-' + df_merged['Hafta'].astype(str).str.zfill(2)
        # Week number to date (Monday of that week)
        df_merged['WeekStartDate'] = df_merged.apply(lambda row: datetime.strptime(f"{row['Yıl']}-W{int(row['Hafta'])}-1", "%Y-W%W-%w").date() if row['Hafta'] != 53 else datetime(row['Yıl'], 12, 28).date() , axis=1) # 53. hafta için basit bir yaklaşım
        # df_merged['WeekStartDate'] = df_merged.apply(lambda row: datetime.strptime(f"{row['Yıl']}-W{int(row['Hafta'])}-1", "%Y-W%W-%w").date(), axis=1) # Bu satır doğru çalışmaz, 53. hafta sorunu olabilir.
        # Yeni tarih hesaplama mantığı (daha sağlam)
        df_merged['WeekStartDate'] = df_merged.apply(lambda row: datetime(row['Yıl'], 1, 1) + timedelta(weeks=row['Hafta']-1), axis=1)
        # Hafta 1'in Pazartesisi'ni bulmak için ilk olarak doğru yılın ilk haftasını bulmalıyız.
        # Week number to date conversion needs to be robust for all years
        # Let's try to reconstruct date from year and week number
        def get_date_from_year_week(year, week):
            try:
                # Assuming week 1 starts on the first Monday of the year
                # This might need adjustment based on ISO week dates vs. arbitrary week starts
                return datetime.strptime(f'{year}-W{int(week)}-1', "%Y-W%W-%w")
            except ValueError:
                # Handle cases like week 53 not existing in certain years, or other parsing errors
                # Fallback to a simpler approach or skip
                return pd.NaT # Not a Time
        
        df_merged['WeekStartDate'] = df_merged.apply(lambda row: get_date_from_year_week(row['Yıl'], row['Hafta']), axis=1)
        df_merged.dropna(subset=['WeekStartDate'], inplace=True) # Geçersiz tarihleri kaldır


        rolling_correlations = []
        for i in range(len(df_merged) - rolling_window_size + 1):
            window_df = df_merged.iloc[i:i + rolling_window_size].copy()
            window_start_date = window_df['WeekStartDate'].min()
            window_end_date = window_df['WeekStartDate'].max()

            # Pencere içinde yeterli veri yoksa atla
            if len(window_df) < rolling_window_size * 0.5: # En az yarısı kadar veri olsun
                continue

            # Her gecikme için korelasyon hesapla
            window_lags_corrs = {'Window_Start_Date': window_start_date, 'Window_End_Date': window_end_date}
            
            for lag in lags:
                if lag == 0:
                    corr, p_val = sps.pearsonr(window_df[virus_column_name], window_df['Haftalık_Sıcaklık'])
                else:
                    df_lagged_window = window_df.copy()
                    df_lagged_window['Lagged_Sıcaklık'] = df_lagged_window['Haftalık_Sıcaklık'].shift(lag)
                    df_lagged_window.dropna(subset=[virus_column_name, 'Lagged_Sıcaklık'], inplace=True)
                    if len(df_lagged_window) > 1:
                        corr, p_val = sps.pearsonr(df_lagged_window[virus_column_name], df_lagged_window['Lagged_Sıcaklık'])
                    else:
                        corr, p_val = np.nan, np.nan
                window_lags_corrs[f'Corr_Lag_{lag}'] = corr
                window_lags_corrs[f'P_Val_Lag_{lag}'] = p_val
            rolling_correlations.append(window_lags_corrs)

        rolling_corr_df = pd.DataFrame(rolling_correlations)
        if not rolling_corr_df.empty:
            print(f"Dönemsel Korelasyon Sonuçları (Pencere Boyutu: {rolling_window_size} Hafta, İlk 5 Satır):\n", rolling_corr_df.head())
            if save_to_file:
                rolling_corr_df.to_excel(os.path.join(output_excel_folder, f'{sanitize_filename(virus_column_name)}_Rolling_Lagged_Correlation_{rolling_window_size}W.xlsx'), index=False)
        else:
            print("UYARI: Dönemsel korelasyon analizi için yeterli veri bulunamadı.")


        # --- Görselleştirmeler (Matplotlib fig objelerini döndüreceğiz) ---

        # 1. Haftalık Prevalans ve Sıcaklık Trendleri (Timeline)
        fig_prevalance_temp_trend, ax = plt.subplots(figsize=(15, 7))
        ax2 = ax.twinx()
        sns.lineplot(x='Date', y=virus_column_name, data=df_merged, ax=ax, label=virus_column_name, color='blue', errorbar=None)
        sns.lineplot(x='Date', y='Haftalık_Sıcaklık', data=df_merged, ax=ax2, label='Haftalık Sıcaklık (°C)', color='red', errorbar=None)
        ax.set_title(f'Haftalık {virus_column_name} ve Sıcaklık Trendleri')
        ax.set_xlabel('Tarih')
        ax.set_ylabel(virus_column_name, color='blue')
        ax2.set_ylabel('Haftalık Sıcaklık (°C)', color='red')
        ax.tick_params(axis='y', labelcolor='blue')
        ax2.tick_params(axis='y', labelcolor='red')
        fig_prevalance_temp_trend.legend(loc='upper left', bbox_to_anchor=(0.1, 1.0))
        plt.tight_layout()
        figures_dict['prevalance_temp_trend'] = fig_prevalance_temp_trend
        if save_to_file:
            fig_prevalance_temp_trend.savefig(os.path.join(output_plots_folder, f'{sanitize_filename(virus_column_name)}_Prevalance_and_Temperature_Trend.png'))
            plt.close(fig_prevalance_temp_trend)


        # 2. Genel Gecikmeli Korelasyon Isı Haritası
        if not overall_corr_df.empty:
            fig_overall_heatmap, ax = plt.subplots(figsize=(8, 6))
            pivot_df = overall_corr_df.pivot_table(index='Gecikme (Hafta)', values='Korelasyon (r)')
            sns.heatmap(pivot_df, annot=True, cmap='coolwarm', fmt=".2f", ax=ax, linewidths=.5, linecolor='black')
            ax.set_title(f'Genel Gecikmeli Korelasyon ({virus_column_name} vs. Haftalık Sıcaklık)')
            plt.tight_layout()
            figures_dict['overall_heatmap'] = fig_overall_heatmap
            if save_to_file:
                fig_overall_heatmap.savefig(os.path.join(output_plots_folder, f'{sanitize_filename(virus_column_name)}_Overall_Lagged_Correlation_Heatmap.png'))
                plt.close(fig_overall_heatmap)


        # 3. Yıllık Gecikmeli Korelasyon Isı Haritası
        if not yearly_corr_df.empty:
            fig_yearly_heatmap, ax = plt.subplots(figsize=(12, len(yearly_corr_df['Yıl'].unique()) * 0.8)) # Yıl sayısına göre boyut ayarla
            yearly_pivot_df = yearly_corr_df.pivot_table(index='Yıl', columns='Gecikme (Hafta)', values='Korelasyon (r)')
            sns.heatmap(yearly_pivot_df, annot=True, cmap='coolwarm', fmt=".2f", ax=ax, linewidths=.5, linecolor='black')
            ax.set_title(f'Yıllara Göre Gecikmeli Korelasyon ({virus_column_name} vs. Haftalık Sıcaklık)')
            plt.tight_layout()
            figures_dict['yearly_heatmap'] = fig_yearly_heatmap
            if save_to_file:
                fig_yearly_heatmap.savefig(os.path.join(output_plots_folder, f'{sanitize_filename(virus_column_name)}_Yearly_Lagged_Correlation_Heatmap.png'))
                plt.close(fig_yearly_heatmap)

        # 4. Dönemsel Korelasyon Trendi (Belirli bir gecikme için)
        if not rolling_corr_df.empty:
            # Örneğin, en güçlü korelasyonun gözlemlendiği gecikme için trend çizelim (burada 0. gecikmeyi seçiyorum, kullanıcıya seçtirilebilir)
            fig_rolling_corr_trend, ax = plt.subplots(figsize=(15, 7))
            
            # Tüm 'Corr_Lag_X' sütunlarını topla
            corr_cols = [col for col in rolling_corr_df.columns if col.startswith('Corr_Lag_')]
            
            if corr_cols: # Eğer korelasyon sütunları varsa
                # Genellikle en ilginç gecikmelerden birini seçeriz, burada 0'ı varsayalım veya kullanıcı seçsin
                lag_to_plot = f'Corr_Lag_0'
                if lag_to_plot in rolling_corr_df.columns:
                    sns.lineplot(x='Window_Start_Date', y=lag_to_plot, data=rolling_corr_df, ax=ax, errorbar=None)
                    ax.set_title(f'{virus_column_name} ve Sıcaklık Arası Dönemsel Korelasyon (Gecikme: 0 Hafta, Pencere: {rolling_window_size} Hafta)')
                    ax.set_xlabel('Pencere Başlangıç Tarihi')
                    ax.set_ylabel('Korelasyon (r)')
                    ax.axhline(0, color='gray', linestyle='--', linewidth=0.8) # 0 çizgisi
                    plt.tight_layout()
                    figures_dict['rolling_corr_trend'] = fig_rolling_corr_trend
                    if save_to_file:
                        fig_rolling_corr_trend.savefig(os.path.join(output_plots_folder, f'{sanitize_filename(virus_column_name)}_Rolling_Correlation_Trend_{rolling_window_size}W.png'))
                        plt.close(fig_rolling_corr_trend)
                else:
                    print(f"UYARI: '{lag_to_plot}' sütunu dönemsel korelasyon DataFrame'inde bulunamadı. Dönemsel trend grafiği çizilemedi.")
            else:
                print("UYARI: Dönemsel korelasyon sütunları bulunamadı. Dönemsel trend grafiği çizilemedi.")

        # 5. Aylık Prevalans Dağılımları (Boxplot)
        df_merged['Ay'] = df_merged['Date'].dt.month
        fig_monthly_boxplot, ax = plt.subplots(figsize=(12, 6))
        sns.boxplot(x='Ay', y=virus_column_name, data=df_merged, ax=ax)
        ax.set_title(f'Aylık {virus_column_name} Dağılımı')
        ax.set_xlabel('Ay')
        ax.set_ylabel(virus_column_name)
        ax.set_xticklabels([calendar.month_abbr[i] for i in range(1, 13)])
        plt.tight_layout()
        figures_dict['monthly_boxplot'] = fig_monthly_boxplot
        if save_to_file:
            fig_monthly_boxplot.savefig(os.path.join(output_plots_folder, f'{sanitize_filename(virus_column_name)}_Monthly_Prevalance_Distribution.png'))
            plt.close(fig_monthly_boxplot)


        # 6. Prevalans vs. Sıcaklık Serpme Grafiği
        fig_scatter_plot, ax = plt.subplots(figsize=(10, 8))
        sns.scatterplot(x='Haftalık_Sıcaklık', y=virus_column_name, data=df_merged, ax=ax, alpha=0.6)
        sns.regplot(x='Haftalık_Sıcaklık', y=virus_column_name, data=df_merged, ax=ax, scatter=False, color='red') # Regresyon çizgisi
        ax.set_title(f'{virus_column_name} vs. Haftalık Sıcaklık Serpme Grafiği')
        ax.set_xlabel('Haftalık Sıcaklık (°C)')
        ax.set_ylabel(virus_column_name)
        plt.tight_layout()
        figures_dict['scatter_plot'] = fig_scatter_plot
        if save_to_file:
            fig_scatter_plot.savefig(os.path.join(output_plots_folder, f'{sanitize_filename(virus_column_name)}_Prevalance_vs_Temperature_Scatter.png'))
            plt.close(fig_scatter_plot)

        if save_to_file:
            print(f"SUCCESS: Tüm görselleştirmeler başarıyla '{output_plots_folder}' klasörüne kaydedildi.")

        return (overall_corr_df, yearly_corr_df, rolling_corr_df, figures_dict)

    except Exception as e:
        print(f"HATA: Gecikmeli korelasyon analizi sırasında genel bir hata oluştu: {e}")
        import traceback
        traceback.print_exc() # Hata izini yazdır
        return (None, None, None, None)

# __name__ == "__main__": bloğunu olduğu gibi bırakın, içindeki test çağrıları Streamlit'te çalışmayacaktır.
# Ancak, fonksiyonumuz artık daha esnek olduğu için bu blok yine de işlevini korur.