
from datetime import datetime
from IPython.display import display, HTML
from plotly import express as px
from plotly import graph_objects as go
import pandas as pd
import numpy as np

import os
# os.environ["PATH"] = "/home/michalula/.cargo/bin:" + os.environ["PATH"]
# ! which modkit
# ! modkit --version
# /home/michalula/.cargo/bin/modkit
# modkit 0.5.1

# ! python3 -m pip install plotly
# ! python3 -m pip install matplotlib
# ! python3 -m pip install nbformat>=4.2.0

def current_time():
    """Returns the current date and time as a formatted string."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S") 

# print("Current Date and Time:", current_time())import os


def load_pileup_bed(bed_path):
    # bed_path = existing[0]
    print("Reading bedMethyl file:", bed_path)

    # bedMethyl column names (18 columns as provided)
    colnames = [
        "chrom", "start", "end", "mod_code", "score", "strand",
        "start2", "end2", "color",
        "Nvalid_cov", "percent_modified", "Nmod", "Ncanonical",
        "Nother_mod", "Ndelete", "Nfail", "Ndiff", "Nnocall"
    ]

    # Configure dtypes where reasonable
    dtypes = {
        "chrom": str,
        "start": "Int64",
        "end": "Int64",
        "mod_code": str,
        "score": "Int64",
        "strand": str,
        "start2": "Int64",
        "end2": "Int64",
        "color": str,
        "Nvalid_cov": "Int64",
        "percent_modified": float,
        "Nmod": "Int64",
        "Ncanonical": "Int64",
        "Nother_mod": "Int64",
        "Ndelete": "Int64",
        "Nfail": "Int64",
        "Ndiff": "Int64",
        "Nnocall": "Int64"
    }

    compression = "gzip" if bed_path.endswith(".gz") else None

    # Read file (headerless BED-like table). If file has extra columns, keep them with automatic numeric conversion below.
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        comment="#",
        names=colnames,
        dtype=dtypes,
        compression=compression,
        engine="python",
        na_values=[".", "NA", ""],
        keep_default_na=True
    )

    # If file contained more than 18 columns, pandas assigned remaining data to extra columns named like col_18, col_19...
    # Ensure numeric conversion for numeric-like columns
    for c in df.columns:
        if df[c].dtype == object:
            # try safe numeric conversion where appropriate
            try:
                df[c] = pd.to_numeric(df[c], errors="ignore")
            except Exception:
                pass

    print("Loaded DataFrame shape:", df.shape)
    display(df.head())
    return df


def plot_pileup_df(df_roi, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    # ensure numeric types for plotting
    df_roi['pos'] = df_roi['start'].astype(int)
    df_roi['percent_modified'] = df_roi['percent_modified'].astype(float)
    df_roi['Nvalid_cov'] = df_roi['Nvalid_cov'].astype(int)
    df_roi['Nmod'] = df_roi['Nmod'].astype(int)
    df_roi['Ncanonical'] = df_roi['Ncanonical'].astype(int)

    # Scatter: genomic position vs percent modified (point size = coverage)
    fig1 = px.scatter(
        df_roi,
        x='pos',
        y='percent_modified',
        color='strand',
        size='Nvalid_cov',
        hover_data=['Nvalid_cov','Nmod','Ncanonical','Nother_mod','Nnocall'],
        title='Percent modified across ROI (size = Nvalid_cov)',
        height=500
    )
    fig1.update_layout(xaxis_title='Genomic position (start)', yaxis_title='Percent modified')
    fig1.show()
    # fig1.write_html(os.path.join(out_dir, "roi_percent_modified_scatter.html"), include_plotlyjs='cdn')

    # Histogram: coverage distribution
    fig2 = px.histogram(
        df_roi,
        x='Nvalid_cov',
        nbins=40,
        title='Distribution of Nvalid_cov (coverage) in ROI',
        height=400
    )
    fig2.update_layout(xaxis_title='Nvalid_cov', yaxis_title='Count')
    fig2.show()
    # fig2.write_html(os.path.join(out_dir, "roi_nvalidcov_hist.html"), include_plotlyjs='cdn')

    # Bar: top sites by percent_modified showing Nmod vs Ncanonical (stacked)
    topn = 30
    df_top = df_roi.sort_values('percent_modified', ascending=False).head(topn).copy()
    df_top = df_top.assign(label=df_top['pos'].astype(str) + ":" + df_top['strand'])
    fig3 = go.Figure()
    fig3.add_trace(go.Bar(name='Nmod', x=df_top['label'], y=df_top['Nmod']))
    fig3.add_trace(go.Bar(name='Ncanonical', x=df_top['label'], y=df_top['Ncanonical']))
    fig3.update_layout(barmode='stack', title=f'Sorted Top {topn} sites by percent_modified (stacked Nmod / Ncanonical)',
                    xaxis_title='position:strand', yaxis_title='reads', height=520)
    fig3.show()
    # fig3.write_html(os.path.join(out_dir, "roi_top_sites_stacked_counts.html"), include_plotlyjs='cdn')

    # Print simple summaries
    print("# Rows:", df_roi.shape[0])
    print("Percent modified: median={:.2f}, mean={:.2f}".format(df_roi['percent_modified'].median(), df_roi['percent_modified'].mean()))
    print("Coverage (Nvalid_cov): min={}, median={}, max={}".format(df_roi['Nvalid_cov'].min(), df_roi['Nvalid_cov'].median(), df_roi['Nvalid_cov'].max()))

    # Display first rows table for quick inspection
    display(HTML(df_roi.head(20).to_html(index=False)))


    # Bar: top sites by percent_modified showing Nmod vs Ncanonical (stacked)
    topn = df_roi.shape[0]
    # df_top = df_roi.sort_values('percent_modified', ascending=False).head(topn).copy()
    df_top = df_roi.copy() #.sort_values('percent_modified', ascending=False).head(topn).copy()
    df_top = df_top.assign(label=df_top['pos'].astype(str) + ":" + df_top['strand'])
    fig3 = go.Figure()
    fig3.add_trace(go.Bar(name='Nmod', x=df_top['label'], y=df_top['Nmod']))
    fig3.add_trace(go.Bar(name='Ncanonical', x=df_top['label'], y=df_top['Ncanonical']))
    fig3.update_layout(barmode='stack', title=f'All {topn} CpG sites by percent_modified (stacked Nmod / Ncanonical) [ordered=not s]',
                    xaxis_title='position:strand', yaxis_title='reads', height=520)
    fig3.show()
    # fig3.write_html(os.path.join(out_dir, "roi_top_sites_stacked_counts.html"), include_plotlyjs='cdn')

    # Bar: top sites by percent_modified showing Nmod vs Ncanonical (stacked, NOT SORTED)
    topn = df_roi.shape[0]
    df_top = df_roi.sort_values('percent_modified', ascending=False).head(topn).copy()
    df_top = df_top.assign(label=df_top['pos'].astype(str) + ":" + df_top['strand'])
    fig3 = go.Figure()
    fig3.add_trace(go.Bar(name='Nmod', x=df_top['label'], y=df_top['Nmod']))
    fig3.add_trace(go.Bar(name='Ncanonical', x=df_top['label'], y=df_top['Ncanonical']))
    fig3.update_layout(barmode='stack', title=f'Sorted Top {topn} sites by percent_modified (stacked Nmod / Ncanonical)',
                    xaxis_title='position:strand', yaxis_title='reads', height=520)
    fig3.show()
    # fig3.write_html(os.path.join(out_dir, "roi_top_sites_stacked_counts.html"), include_plotlyjs='cdn')

    # Print simple summaries
    print("# rows:", df_roi.shape[0])
    print("Percent modified: median={:.2f}, mean={:.2f}".format(df_roi['percent_modified'].median(), df_roi['percent_modified'].mean()))
    print("Coverage (Nvalid_cov): min={}, median={}, max={}".format(df_roi['Nvalid_cov'].min(), df_roi['Nvalid_cov'].median(), df_roi['Nvalid_cov'].max()))

    # Display first rows table for quick inspection
    display(HTML(df_roi.head(20).to_html(index=False)))

    # Bar: top sites by percent_modified showing Nmod vs Ncanonical (stacked) percentages
    topn = 30
    df_top = df_roi.sort_values('percent_modified', ascending=False).head(topn).copy()
    df_top = df_top.assign(label=df_top['pos'].astype(str) + ":" + df_top['strand'])
    df_top = df_top.assign(Ntotal=df_top['Nmod'] + df_top['Ncanonical'])
    df_top = df_top.assign(Nmod_perc=(df_top['Nmod'] / df_top['Ntotal']) * 100)
    df_top = df_top.assign(Ncanonical_perc=(df_top['Ncanonical'] / df_top['Ntotal']) * 100)
    fig4 = go.Figure()
    fig4.add_trace(go.Bar(name='Nmod %', x=df_top['label'], y=df_top['Nmod_perc']))
    fig4.add_trace(go.Bar(name='Ncanonical %', x=df_top['label'], y=df_top['Ncanonical_perc']))
    fig4.update_layout(barmode='stack', title=f'Top {topn} sites by percent_modified (stacked Nmod % / Ncanonical %)',
                    xaxis_title='position:strand', yaxis_title='percentage', height=520)
    fig4.show()
    # fig4.write_html(os.path.join(out_dir, "roi_top_sites_stacked_percentage.html"), include_plotlyjs='cdn')     


    # Bar: top sites by percent_modified showing Nmod vs Ncanonical (stacked) percentages
    topn = 277
    df_top = df_roi.sort_values('percent_modified', ascending=False).head(topn).copy()
    df_top = df_top.assign(label=df_top['pos'].astype(str) + ":" + df_top['strand'])
    df_top = df_top.assign(Ntotal=df_top['Nmod'] + df_top['Ncanonical'])
    df_top = df_top.assign(Nmod_perc=(df_top['Nmod'] / df_top['Ntotal']) * 100)
    df_top = df_top.assign(Ncanonical_perc=(df_top['Ncanonical'] / df_top['Ntotal']) * 100)
    fig4 = go.Figure()
    fig4.add_trace(go.Bar(name='Nmod %', x=df_top['label'], y=df_top['Nmod_perc']))
    fig4.add_trace(go.Bar(name='Ncanonical %', x=df_top['label'], y=df_top['Ncanonical_perc']))
    fig4.update_layout(barmode='stack', title=f'Top {topn} sites by percent_modified (stacked Nmod % / Ncanonical %)',
                    xaxis_title='position:strand', yaxis_title='percentage', height=520)
    fig4.show()
    # fig4.write_html(os.path.join(out_dir, "roi_top_sites_stacked_percentage.html"), include_plotlyjs='cdn')     

    # Bar: Unsorted sites by percent_modified showing Nmod vs Ncanonical (stacked) percentages
    df_top = df_roi.copy()
    df_top = df_top.assign(label=df_top['pos'].astype(str) + ":" + df_top['strand'])
    df_top = df_top.assign(Ntotal=df_top['Nmod'] + df_top['Ncanonical'])
    df_top = df_top.assign(Nmod_perc=(df_top['Nmod'] / df_top['Ntotal']) * 100)
    df_top = df_top.assign(Ncanonical_perc=(df_top['Ncanonical'] / df_top['Ntotal']) * 100)
    fig5 = go.Figure()
    fig5.add_trace(go.Bar(name='Nmod %', x=df_top['label'], y=df_top['Nmod_perc']))
    fig5.add_trace(go.Bar(name='Ncanonical %', x=df_top['label'], y=df_top['Ncanonical_perc']))
    fig5.update_layout(barmode='stack', title=f'All sites by percent_modified (stacked Nmod % / Ncanonical %)',
                    xaxis_title='position:strand', yaxis_title='percentage', height=520)
    fig5.show()
    # fig5.write_html(os.path.join(out_dir, "roi_all_sites_stacked_percentage.html"), include_plotlyjs='cdn')    


    return df_top


def load_dmr_and_parse(dmr_path, dmr_output_path, date_today):
    out_dir = dmr_output_path
    print("out_dir:", out_dir)

    # Read DMR BED (robust to header/no-header) and assign canonical column names (uses existing vars: dmr_path, out_dir, date_today, pd, os)
    canonical_cols = [
        "chrom", "start", "end", "name", "score", "strand",
        "samplea_counts", "samplea_total", "sampleb_counts", "sampleb_total",
        "samplea_percents", "sampleb_percents",
        "samplea_fraction_modified", "sampleb_fraction_modified",
        "map_pvalue", "effect_size",
        "cohen_h", "cohen_h_low", "cohen_h_high",
    ]
        # "balanced_map_pvalue", "balanced_effect_size"

    # read file with header and fallback to header=None when headers look numeric or columns are unexpected
    try:
        dmr_df = pd.read_csv(dmr_path, sep="\t", comment="#", header=None, engine="python")

        # dmr_df = pd.read_csv(dmr_path, sep="\t", comment="#", engine="python") # , header=0
        # # heuristic: if too many numeric-looking column names, re-read as headerless
        # numeric_headers = sum(1 for c in dmr_df.columns if str(c).strip().isdigit())
        # if numeric_headers >= (len(dmr_df.columns) / 2) or dmr_df.shape[1] < 3:
        #     dmr_df = pd.read_csv(dmr_path, sep="\t", comment="#", header=None, engine="python")
    except Exception:
        dmr_df = pd.read_csv(dmr_path, sep="\t", comment="#", header=None, engine="python")

    # assign canonical names up to number of columns present, add generic names for extras
    ncols = dmr_df.shape[1]
    if ncols <= len(canonical_cols):
        dmr_df.columns = canonical_cols[:ncols]
    else:
        extras = [f"col_{i}" for i in range(ncols - len(canonical_cols))]
        dmr_df.columns = canonical_cols + extras

    # coerce obvious numeric columns to numeric where present
    num_cols_to_try = [
        "start", "end", "score",
        "samplea_total", "sampleb_total",
        "samplea_fraction_modified", "sampleb_fraction_modified",
        "map_pvalue", "effect_size",
        "balanced_map_pvalue", "balanced_effect_size"
    ]
    for c in num_cols_to_try:
        if c in dmr_df.columns:
            dmr_df[c] = pd.to_numeric(dmr_df[c], errors="coerce")

    # ensure output directory exists and save parsed table (parquet preferred)
    os.makedirs(out_dir, exist_ok=True)
    parsed_path = os.path.join(out_dir, f"{date_today}_dmr_parsed.parquet")
    try:
        dmr_df.to_parquet(parsed_path, index=False)
        print("Saved parquet:", parsed_path)
    except Exception:
        csv_path = os.path.join(out_dir, f"{date_today}_dmr_parsed.csv")
        dmr_df.to_csv(csv_path, index=False)
        print("Parquet not available, saved csv:", csv_path)

    print("Loaded DMR:", dmr_path)
    print("Assigned columns:", dmr_df.columns.tolist())
    print("Shape:", dmr_df.shape)
    return dmr_df


def plot_dmr_summary(dmr_df, dmr_output_path, date_today):
    os.makedirs(dmr_output_path, exist_ok=True)

    # Save a table summary
    summary = dmr_df.describe(include='all').transpose()
    summary_path = os.path.join(dmr_output_path, f"{date_today}_dmr_column_summary.csv")
    summary.to_csv(summary_path)

    numcols = dmr_df.select_dtypes(include=['number']).columns.tolist()

    def _safe_name(name):
        return str(name).replace(os.sep, "_").replace(" ", "_").replace("\t", "_")

    # Per-column visualizations
    for col in dmr_df.columns:
        safe = _safe_name(col)
        try:
            if col in numcols:
                # Histogram
                fig_h = px.histogram(dmr_df, x=col, nbins=80, title=f"Histogram: {col}")
                # fig_h.write_html(os.path.join(out_dir, f"{date_today}_dmr_hist_{safe}.html"), include_plotlyjs='cdn')
                fig_h.show()

                # Boxplot
                fig_b = px.box(dmr_df, y=col, points="outliers", title=f"Boxplot: {col}")
                # fig_b.write_html(os.path.join(out_dir, f"{date_today}_dmr_box_{safe}.html"), include_plotlyjs='cdn')
                fig_b.show()
            else:
                # Categorical / text: show top value counts (up to 50)
                vc = dmr_df[col].fillna("NA").astype(str).value_counts().head(50)
                if len(vc):
                    fig_c = px.bar(x=vc.values[::-1], y=vc.index.astype(str)[::-1], orientation='h',
                                title=f"Top value counts: {col}", labels={'x':'count','y':col})
                    fig_c.update_layout(yaxis={'categoryorder':'array','categoryarray':vc.index[::-1].astype(str).tolist()})
                    # fig_c.write_html(os.path.join(out_dir, f"{date_today}_dmr_valcounts_{safe}.html"), include_plotlyjs='cdn')
                    fig_c.show()
                else:
                    # fallback: display empty info
                    display(HTML(f"<b>{col}</b>: no values to plot"))
        except Exception as e:
            print(f"Skipped plotting column {col!r} due to error: {e}")

    # Correlation heatmap for numeric columns
    if len(numcols) >= 2:
        try:
            corr = dmr_df[numcols].corr()
            fig_corr = px.imshow(corr, text_auto=True, aspect="auto", title="Correlation matrix (numeric columns)")
            # fig_corr.write_html(os.path.join(out_dir, f"{date_today}_dmr_correlation_numeric.html"), include_plotlyjs='cdn')
            fig_corr.show()
        except Exception as e:
            print("Failed to create correlation heatmap:", e)

    print("Saved summary:", summary_path)
    print("Plots saved to:", dmr_output_path)

# GET ROI data subset:
    # 137*2, 277-5
    # pileup_K562_Zoff_Filter_pileup_df[pileup_K562_Zoff_Filter_pileup_df['start'] == 206583387]
    # pileup_K562_Zoff_Filter_pileup_df[pileup_K562_Zoff_Filter_pileup_df['start'] == 206583388]
    # pileup_K562_Zoff_Filter_pileup_df[pileup_K562_Zoff_Filter_pileup_df['start'] == 206589746]
    # (279-69) / 2
    # pileup_K562_Zoff_Filter_pileup_df_roi = pileup_K562_Zoff_Filter_pileup_df.iloc[69:343, :]  # Display target region rows
    # print(pileup_K562_Zoff_Filter_pileup_df_roi.shape, pileup_K562_Zoff_Filter_pileup_df_roi.shape[0]/2)
    # pileup_K562_Zoff_Filter_pileup_df_roi