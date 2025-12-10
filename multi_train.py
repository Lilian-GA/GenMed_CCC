import os
import scanpy as sc
import cell2location as c2l
import multiprocessing as mp
import commot as ct
import squidpy as sq
def run_cell2location(sc_path, st_path, outdir="models"):
    # --- Extraction d’un nom propre basé sur le fichier ---
    sc_name = os.path.basename(sc_path).replace(".h5ad", "")
    print(sc_name)
    st_name = os.path.basename(st_path).replace(".h5ad", "")
    base_name = f"{sc_name}__{st_name}"
    outdata_dir = os.path.join("data","data_processed")
    os.makedirs(outdir, exist_ok=True)

    ######
    print("########################################################################")
    print("########################################################################")


    # --- Lecture des données ---
    adata_sc = sc.read(sc_path)
    adata_st = sc.read(st_path)

    # Spatial preprocessing
    adata_st.obsm['spatial'] = adata_st.obsm.get("X_spatial", adata_st.obsm.get("spatial"))
    adata_st.layers["counts"] = adata_st.X.copy()
    print(f"[INFO] Traitement du couple : {base_name}")


    adata_st.var["MT_gene"] = [g.startswith("MT-") for g in adata_st.var["features"]]
    adata_st.obsm["MT"] = adata_st[:, adata_st.var["MT_gene"].values].X.toarray()
    adata_st = adata_st[:, ~adata_st.var["MT_gene"].values]

    # --- Shared features ---
    shared = [g for g in adata_st.var_names if g in adata_sc.var_names]
    adata_sc = adata_sc[:, shared].copy()
    adata_st = adata_st[:, shared].copy()
    print("------------ adata_sc ----------")
    print(adata_sc)
    print("------------ adata_st ----------")
    print(adata_st)
    #
    # # --- Modèle scRNAseq ---
    # c2l.models.RegressionModel.setup_anndata(
    #     adata=adata_sc,
    #     batch_key="Piece_ID",
    #     labels_key="cell_type",
    # )
    #
    # model = c2l.models.RegressionModel(adata_sc)
    # print("--------------- Training single cell reference -----------------------")
    # model.train(max_epochs=750, batch_size=2500, train_size=1, lr=0.01)
    #
    # model_path = os.path.join(outdir, f"{base_name}__sc_regression")
    # model.save(model_path, overwrite=True)
    #
    # model.export_posterior(
    #     adata_sc,
    #     sample_kwargs={"num_samples": 1000, "batch_size": 2500},
    # )
    #
    # inf_aver = adata_sc.varm['means_per_cluster_mu_fg']
    #
    # # --- Modèle spatial ---
    # print("--------------- deconvolution -----------------------")
    # c2l.models.Cell2location.setup_anndata(adata_st, layers = "counts")
    # mod_sp = c2l.models.Cell2location(
    #     adata_st,
    #     cell_state_df=inf_aver,
    #     N_cells_per_location=8
    # )
    # mod_sp.train(max_epochs=800, lr = 0.01)
    #
    # mod_sp_path = os.path.join(outdir, f"{base_name}__spatial_model")
    # mod_sp.save(mod_sp_path)
    #
    # print(f"[DONE] Modèles enregistrés pour {base_name}")
    #
    # # Exporter abondances (q05 = confiance haute)
    # adata_st = mod_sp.export_posterior(adata_st,
    #                                    sample_kwargs={'num_samples': 1000})
    # adata_st.obs[mod_sp.factor_names_] = adata_st.obsm['q05_cell_abundance_w_sf']
    #
    # # Clustering régions tissulaires
    # sc.pp.neighbors(adata_st, use_rep='q05_cell_abundance_w_sf')
    # sc.tl.leiden(adata_st)
    # sc.tl.umap(adata_st)
    # sq.gr.spatial_neighbors(adata_st)
    #
    # # cell-cell communication
    # ct.tl.spatial_communication(adata_st,
    #                             database_name=base_name,
    #                             df_ligrec=CELLCHAT, dis_thr=500, heteromeric=True)
    #
    # # 1. Extraire counts bruts
    # raw_counts = adata_st.layers["counts"].copy()
    #
    # # 2. Nettoyer var : retirer `_index` si présent
    # var = adata_st.var.copy()
    # if "_index" in var.columns:
    #     var = var.drop(columns=["_index"])
    #
    # # 3. Créer un nouvel AnnData pour raw
    # adata_st.raw = sc.AnnData(
    #     X=raw_counts,
    #     obs=adata_st.obs.copy(),
    #     var=var
    # )
    #
    # # 1. Extraire counts bruts
    # raw_counts = adata_sc.layers["counts"].copy() if "counts" in adata_sc.layers else adata_sc.X.copy()
    #
    # # 2. Nettoyer var : retirer `_index` si présent
    # var = adata_sc.var.copy()
    # if "_index" in var.columns:
    #     var = var.drop(columns=["_index"])
    #
    # # 3. Recréer un raw propre
    # adata_sc.raw = sc.AnnData(
    #     X=raw_counts,
    #     obs=adata_sc.obs.copy(),
    #     var=var
    # )
    #
    # #Save data
    # data_st_path = os.path.join(outdir, f"{st_name}_st.h5ad")
    # data_sc_path = os.path.join(outdir, f"{sc_name}_sc.h5ad")
    # adata_st.write(data_st_path)
    # adata_sc.write(data_sc_path)
def run_all_pairs(pairs, n_jobs=2):
    with mp.Pool(processes=n_jobs) as pool:
        pool.starmap(run_cell2location, pairs)


if __name__ == "__main__":
    CELLCHAT = ct.pp.ligand_receptor_database(species='human',
                                              database='CellChat')
    # Liste des couples (scRNAseq_file, spatial_file)
    pairs = [
        ("data/data_h5ad/HT224P1-S1.h5ad",
         "data/data_h5ad/HT224P1-S1Fc2U1Z1Bs1-SeuratObj.h5ad"),

        ("data/data_h5ad/HT259P1-S1H1.h5ad",
        "data/data_h5ad/HT259P1-S1H1Fc2U1Z1Bs1-SeuratObj.h5ad"),

        ("data/data_h5ad/snRNA_L4__HT243B1-H3A2.h5ad",
        "data/data_h5ad/HT243B1H4A2-S1Fc1U2Z1B1-SeuratObj.h5ad"),

        ("data/data_h5ad/snRNA_L4__HT425B1-S1H1_combo.h5ad",
        "data/data_h5ad/HT425B1-S1H2Fs1U1Bp1-SeuratObj.h5ad")
    ]
    print("-------------------- Beginning process ------------------")
    run_all_pairs(pairs, n_jobs=2)
    print("---------------- Done process -------------------")
