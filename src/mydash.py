import os
import scanpy as sc
import plotly.express as px
from dash import Dash, dcc, html, Input, Output
import numpy as np

# Load data once, static file name
adata = sc.read_h5ad("zebra_light.h5ad")
# Preprocess gene name lookup dictionary (case-insensitive)
gene_lookup = {g.lower().strip(): g for g in adata.var_names}

# Robust gene finding function
def find_gene(gene_name, lookup_dict):
    if not gene_name:
        return None
    return lookup_dict.get(gene_name.strip().lower(), None)

app = Dash(__name__)

app.layout = html.Div([
    html.H2("Gene expression on UMAP"),
    dcc.Input(id='gene-input', type='text', placeholder='Enter gene name...', style={'width': '300px'}),
    html.Div(id='output-message', style={'marginTop': '10px', 'color': 'red'}),
    dcc.Graph(id='gene-plot')
])

@app.callback(
    Output('gene-plot', 'figure'),
    Output('output-message', 'children'),
    Input('gene-input', 'value')
)
def update_plot(gene_input):
    if not gene_input or gene_input.strip() == "":
        return {}, ""

    found_gene = find_gene(gene_input, gene_lookup)
    if not found_gene:
        return {}, f"Gene '{gene_input}' not found in dataset."

    # Extract gene expression
    expr = adata[:, found_gene].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray().flatten()
    else:
        expr = np.array(expr).flatten()

    fig = px.scatter(
        x=adata.obsm['X_umap'][:, 0],
        y=adata.obsm['X_umap'][:, 1],
        color=expr,
        color_continuous_scale='Viridis',
        labels={'color': found_gene},
        title=f"Expression of {found_gene} on UMAP"
    )

    return fig, ""

if __name__ == '__main__':
    # Get port from environment variable PORT, default to 8050 locally
    port = int(os.environ.get('PORT', 8050))
    app.run(host='0.0.0.0', port=port)

