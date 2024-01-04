PRESETS = {"ont": "map-ont", "illumina": "sr"}
CLASSIFY_DBS = {
    "standard": RESULTS / "db/kraken/standard",
    "standard-8": RESULTS / "db/kraken/standard-8",
    "mycobacterium": RESULTS / "db/GTDB_genus_Mycobacterium/kraken/db/hash.k2d",
}

def infer_kraken_db(wildcards):
    if wildcards.lib == "default":
        return RESULTS / "dehumanise/kraken/db/k35/l31/db"
    elif wildcards.lib == "HPRC":
        return RESULTS / "db/HPRC/kraken/db"
    else:
        raise NotImplementedError(f"Kraken lib {wildcards.lib} not known")


def infer_minimap2_db(wildcards):
    preset = PRESETS[wildcards.tech]
    if wildcards.db == "clockwork":
        return RESULTS / f"classify/minimap2/db/db.{preset}.mmi"
    elif wildcards.db == "mtbc":
        return RESULTS / f"db/GTDB_genus_Mycobacterium/MTB.{preset}.mmi"
    elif wildcards.db == "mycobacterium":
        return RESULTS / f"classify/minimap2/db/Mycobacterium.rep.{preset}.mmi"
    else:
        raise ValueError(f"Don't recognise db {wildcards.db}")


def infer_minimap2_db_real(wildcards):
    preset = PRESETS[wildcards.tech]
    if wildcards.db == "clockwork":
        return RESULTS / f"db/minimap2/db.{preset}.mmi"
    elif wildcards.db == "mtbc":
        return RESULTS / f"db/GTDB_genus_Mycobacterium/MTB.{preset}.mmi"
    elif wildcards.db == "mycobacterium":
        return RESULTS / f"db/GTDB_genus_Mycobacterium/Mycobacterium.rep.{preset}.mmi"
    else:
        raise ValueError(f"Don't recognise db {wildcards.db}")


def infer_classify_reads(wildcards):
    if wildcards.tech == "ont":
        return RESULTS / "dehumanise/metagenome.dehumanised.ont.fq.gz"
    elif wildcards.tech == "illumina":
        return [
            RESULTS / f"dehumanise/metagenome_R{i}.dehumanised.illumina.fq.gz"
            for i in [1, 2]
        ]
    else:
        raise ValueError(f"Don't recognise tech {wildcards.tech}")


def infer_classify_real_reads(wildcards):
    if wildcards.tech == "ont":
        return RESULTS / "real/reads/nonhuman.ont.fq.gz"
    elif wildcards.tech == "illumina":
        return [RESULTS / f"real/reads/nonhuman_{i}.illumina.fq.gz" for i in [1, 2]]
    else:
        raise ValueError(f"Don't recognise tech {wildcards.tech}")
