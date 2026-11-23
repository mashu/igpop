//! Step 4: the output tables. Every table is declared as a list of columns,
//! each a name with the function that fills it, so a header and its rows
//! cannot drift apart. The meaning of every column is in docs/src/outputs.md.

use crate::assess::{AlleleResult, Class, ScoreModel, VariantResult};
use crate::evidence::{absent_upper, Evidence};
use crate::freq::{self, PopCount, GNOMAD_GROUPS, POPS};
use crate::hap;
use crate::locus::{Called, Placed, Placement, Window};
use crate::{render, score};
use anyhow::Result;
use std::fmt::Write;

type Cells<T> = Box<dyn Fn(&T) -> Vec<String>>;

/// A tab-separated table over rows of type `T`.
pub struct Table<T> {
    names: Vec<String>,
    blocks: Vec<(usize, Cells<T>)>,
}

impl<T> Default for Table<T> {
    fn default() -> Self {
        Table {
            names: Vec::new(),
            blocks: Vec::new(),
        }
    }
}

impl<T> Table<T> {
    pub fn new() -> Self {
        Self::default()
    }

    /// One column.
    fn col(mut self, name: &str, cell: impl Fn(&T) -> String + 'static) -> Self {
        self.names.push(name.to_string());
        self.blocks.push((1, Box::new(move |row| vec![cell(row)])));
        self
    }

    /// Several columns filled together (one source, one case distinction).
    fn block(mut self, names: Vec<String>, cells: impl Fn(&T) -> Vec<String> + 'static) -> Self {
        self.blocks.push((names.len(), Box::new(cells)));
        self.names.extend(names);
        self
    }

    pub fn header(&self) -> String {
        self.names.join("\t")
    }

    pub fn row(&self, row: &T) -> String {
        let mut out: Vec<String> = Vec::with_capacity(self.names.len());
        for (width, cells) in &self.blocks {
            let cells = cells(row);
            assert_eq!(cells.len(), *width, "a column block must fill all its columns");
            out.extend(cells);
        }
        out.join("\t")
    }

    /// The whole table as text: a header line, then one line per row.
    pub fn render<'a>(&self, rows: impl IntoIterator<Item = &'a T>) -> String
    where
        T: 'a,
    {
        let mut out = self.header();
        out.push('\n');
        for row in rows {
            out.push_str(&self.row(row));
            out.push('\n');
        }
        out
    }
}

/// A number: `.` if undefined, scientific notation if very small.
pub fn num(x: f64) -> String {
    if x.is_nan() {
        ".".into()
    } else if x != 0.0 && x.abs() < 1e-4 {
        format!("{x:.3e}")
    } else {
        format!("{x:.6}")
    }
}

fn dots(n: usize) -> Vec<String> {
    vec![".".to_string(); n]
}

fn list(ids: &[String]) -> String {
    if ids.is_empty() {
        ".".to_string()
    } else {
        ids.join(";")
    }
}

fn names(columns: &[&str]) -> Vec<String> {
    columns.iter().map(|c| c.to_string()).collect()
}

fn per_pop(prefix: &str) -> Vec<String> {
    POPS.iter().map(|p| format!("{prefix}{p}")).collect()
}

fn counts(pops: &[PopCount]) -> Vec<String> {
    pops.iter().map(|c| format!("{}/{}", c.ac, c.an)).collect()
}

/// One row of `variants.tsv`.
pub struct VariantRow<'a> {
    pub ev: &'a Evidence,
    pub allele: &'a AlleleResult<'a>,
    pub called: &'a Called,
    pub v: &'a VariantResult,
    /// `--prior`, for the probability column.
    pub prior: Option<f64>,
}

pub fn variants_table<'a>() -> Table<VariantRow<'a>> {
    type R<'a> = VariantRow<'a>;
    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
    Table::new()
        .col("allele_name", |r: &R| r.allele.placed.allele.clone())
        .col("gene", |r: &R| r.allele.placed.gene.clone())
        .col("rsid", |r: &R| r.v.rsid.clone())
        .col("evidence", |r: &R| {
            if r.v.observed.is_empty() {
                "-".to_string()
            } else {
                r.v.observed.join(";")
            }
        })
        .col("contig", |r: &R| r.allele.placed.contig().to_string())
        .col("genome_pos", |r: &R| r.called.key.0.to_string())
        .col("aln_pos", |r: &R| r.called.allele_pos.to_string())
        .col("type", |r: &R| r.called.kind.as_str().to_string())
        .col("query_variant", move |r: &R| text(&r.called.key.2))
        .col("vcf_ref", |r: &R| {
            r.ev.kgp.get(&r.v.lookup).map_or("-".to_string(), |s| s.vcf_ref.clone())
        })
        .col("vcf_alt", |r: &R| {
            r.ev.kgp.get(&r.v.lookup).map_or("-".to_string(), |s| s.vcf_alt.clone())
        })
        .col("verdict", |r: &R| r.v.verdict.as_str().to_string())
        .col("filters", |r: &R| {
            if r.v.filters.is_empty() {
                ".".to_string()
            } else {
                r.v.filters.iter().map(|f| f.as_str()).collect::<Vec<_>>().join(",")
            }
        })
        .col("score", |r: &R| r.v.score.map_or(".".to_string(), |x| format!("{x:.3}")))
        .col("probability_real", |r: &R| match (r.prior, r.v.score) {
            (Some(prior), Some(x)) => format!("{:.3}", score::posterior(prior, x)),
            _ => ".".to_string(),
        })
        .col("support", |r: &R| if r.ev.counted() { num(r.v.support) } else { ".".to_string() })
        .col("class", |r: &R| r.v.class.as_str().to_string())
        .col("class_basis", |r: &R| r.v.basis.clone().unwrap_or_else(|| ".".to_string()))
        .col("class_p", |r: &R| r.v.class_p.map_or(".".to_string(), num))
        .block(kgp_names(), kgp_cells)
        .block(
            names(&[
                "gnomad_filter",
                "gnomad_flags",
                "gnomad_AC",
                "gnomad_AN",
                "gnomad_AF",
                "gnomad_nhomalt",
                "gnomad_max_group",
                "gnomad_max_group_AF",
                "gnomad_max_group_AF_lower95",
                "gnomad_any_group_AF_upper95",
            ]),
            gnomad_cells,
        )
        .block(
            names(&["dbsnp_common", "dbsnp_datasets", "dbsnp_max_dataset", "dbsnp_max_AF", "dbsnp_freq"]),
            dbsnp_cells,
        )
        .col("locus_status", |r: &R| r.allele.placed.status.as_str().to_string())
        .col("shares_locus_with", |r: &R| r.allele.placed.shared())
}

fn kgp_names() -> Vec<String> {
    let mut out = names(&["kgp_site_qc", "kgp_filter", "kgp_AC", "kgp_AN", "kgp_AF"]);
    out.extend(per_pop("kgp_AF_"));
    out.extend(per_pop("kgp_AC/AN_"));
    out.extend(names(&[
        "kgp_max_pop",
        "kgp_max_pop_AF",
        "kgp_max_pop_AF_lower95",
        "kgp_any_pop_AF_upper95",
        "kgp_pop_heterogeneity_p",
        "kgp_het",
        "kgp_hom",
        "kgp_excess_het_p",
        "kgp_mendelian_error_rate",
        "kgp_vqslod",
    ]));
    out
}

/// 1000 Genomes 30x, unrelated samples. Three cases: a record exists; the
/// data was read but has no record (the allele is absent among the
/// chromosomes of a fully called site); the data was not read.
fn kgp_cells(r: &VariantRow) -> Vec<String> {
    match r.ev.kgp.get(&r.v.lookup) {
        Some(site) => {
            let total = site.total();
            let (max_pop, strongest) = site.strongest_pop();
            let mut out = vec![
                if r.v.site_qc.is_empty() {
                    "PASS".to_string()
                } else {
                    r.v.site_qc.join(",")
                },
                site.filter.clone(),
                total.ac.to_string(),
                total.an.to_string(),
                num(total.af()),
            ];
            out.extend(site.pops.iter().map(|c| num(c.af())));
            out.extend(counts(&site.pops));
            out.extend([
                max_pop.to_string(),
                num(strongest.af()),
                num(strongest.lower95()),
                num(site.max_upper95()),
                num(site.heterogeneity_p()),
                total.het.to_string(),
                total.hom.to_string(),
                num(site.exc_het),
                num(site.mendel),
                num(site.vqslod),
            ]);
            out
        }
        None if r.ev.kgp_read() => {
            let mut out = vec![
                ".".to_string(),
                ".".to_string(),
                "0".to_string(),
                r.ev.kgp_full_total().to_string(),
                num(0.0),
            ];
            out.extend(vec![num(0.0); 5]);
            out.extend(r.ev.kgp_full.iter().map(|an| format!("0/{an}")));
            out.extend(dots(3));
            out.push(num(absent_upper(&r.ev.kgp_full)));
            out.extend(dots(6));
            out
        }
        None => dots(25),
    }
}

/// gnomAD genomes; the same three cases.
fn gnomad_cells(r: &VariantRow) -> Vec<String> {
    match r.ev.gnomad.get(&r.v.lookup) {
        Some(site) => {
            let i = freq::strongest(&site.groups);
            let top = site.groups[i];
            vec![
                site.filter.clone(),
                site.flags.clone(),
                site.total.ac.to_string(),
                site.total.an.to_string(),
                num(site.total.af()),
                site.total.hom.to_string(),
                GNOMAD_GROUPS[i].to_string(),
                num(top.af()),
                num(top.lower95()),
                num(freq::max_upper95(&site.groups)),
            ]
        }
        None if r.ev.gnomad_read() => {
            let mut out = vec![
                ".".to_string(),
                ".".to_string(),
                "0".to_string(),
                r.ev.gnomad_full_total.to_string(),
                num(0.0),
                "0".to_string(),
            ];
            out.extend(dots(3));
            out.push(num(absent_upper(&r.ev.gnomad_full)));
            out
        }
        None => dots(10),
    }
}

/// dbSNP: the COMMON flag for exactly this allele, and the frequencies
/// relayed from the datasets in use.
fn dbsnp_cells(r: &VariantRow) -> Vec<String> {
    if r.ev.sources.dbsnp.is_none() {
        return dots(5);
    }
    let common = if r.ev.dbsnp.common.contains(&r.v.lookup) { "yes" } else { "no" };
    let relayed = r.ev.dbsnp.relayed(&r.v.lookup);
    match relayed.first() {
        Some((top, af)) => {
            let all: Vec<String> = relayed.iter().map(|(d, af)| format!("{d}:{}", num(*af))).collect();
            vec![common.to_string(), relayed.len().to_string(), top.clone(), num(*af), all.join(";")]
        }
        None => vec![
            common.to_string(),
            "0".to_string(),
            ".".to_string(),
            ".".to_string(),
            ".".to_string(),
        ],
    }
}

/// `alleles.tsv` and `final.tsv`: one row per allele.
pub fn alleles_table<'a>() -> Table<AlleleResult<'a>> {
    type R<'a> = AlleleResult<'a>;
    let of = |r: &R, class: Class| r.ids(|v| v.class == class);
    // "Uncommon" keeps IgSNPer's meaning: every variant not classed common.
    let not_common = move |r: &R| -> Vec<String> {
        [Class::Uncommon, Class::Absent, Class::NotAssessed]
            .iter()
            .flat_map(|c| of(r, *c))
            .collect()
    };
    let missing = |r: &R, end: usize| match &r.placed.placement {
        Some(p) => [p.missing.0, p.missing.1][end].to_string(),
        None => ".".to_string(),
    };
    Table::new()
        .col("allele_name", |r: &R| r.placed.allele.clone())
        .col("gene", |r: &R| r.placed.gene.clone())
        .col("verdict", |r: &R| r.verdict.as_str().to_string())
        .col("filters", |r: &R| r.filters())
        .col("reason", |r: &R| r.reason.clone())
        .col(
            "notes",
            |r: &R| if r.notes.is_empty() { ".".to_string() } else { r.notes.join("; ") },
        )
        .col("hap_status", |r: &R| r.hap_status().to_string())
        .col("locus_status", |r: &R| r.placed.status.as_str().to_string())
        .col("shares_locus_with", |r: &R| r.placed.shared())
        .col("trusted", |r: &R| if r.calibrates() { "yes" } else { "no" }.to_string())
        .col("missing_5prime", move |r: &R| missing(r, 0))
        .col("missing_3prime", move |r: &R| missing(r, 1))
        .col("n_variants", |r: &R| r.variants.len().to_string())
        .col("min_score", |r: &R| r.weakest().map_or(".".to_string(), |(x, _)| format!("{x:.3}")))
        .col("min_score_variant", |r: &R| {
            r.weakest().map_or(".".to_string(), |(_, id)| id.to_string())
        })
        .col("common", move |r: &R| of(r, Class::Common).len().to_string())
        .col("common_str", move |r: &R| list(&of(r, Class::Common)))
        .col("uncommon", move |r: &R| not_common(r).len().to_string())
        .col("uncommon_str", move |r: &R| list(&not_common(r)))
        .col("absent", move |r: &R| of(r, Class::Absent).len().to_string())
        .col("absent_str", move |r: &R| list(&of(r, Class::Absent)))
        .col("not_assessed", move |r: &R| of(r, Class::NotAssessed).len().to_string())
        .col("site_qc_flagged", |r: &R| r.ids(|v| !v.site_qc.is_empty()).len().to_string())
        .col("site_qc_flagged_str", |r: &R| list(&r.ids(|v| !v.site_qc.is_empty())))
        .block(haplotype_names(), haplotype_cells)
}

fn haplotype_names() -> Vec<String> {
    let mut out = names(&[
        "hap_tested",
        "hap_untested_str",
        "hap_total",
        "hap_exact",
        "hap_exact_AF",
        "hap_exact_AF_lower95",
        "hap_exact_AF_upper95",
    ]);
    out.extend(per_pop("hap_exact_"));
    out.extend(names(&[
        "hap_exact_max_pop",
        "hap_exact_pop_heterogeneity_p",
        "hap_carrying",
        "hap_carrying_AF",
    ]));
    out.extend(per_pop("hap_carrying_"));
    out
}

/// The haplotype check. `hap_exact_*` only when every variant was tested;
/// `hap_carrying_*` only when at least one was (or the allele has none).
fn haplotype_cells(r: &AlleleResult) -> Vec<String> {
    let Some(h) = &r.haplotype else {
        return dots(haplotype_names().len());
    };
    let total = r.variants.len();
    let untested: Vec<String> = h.untested.iter().map(|i| r.variants[*i].id.clone()).collect();
    let founders = hap::total(&h.carrying).an;
    let mut out = vec![format!("{}/{total}", total - untested.len()), list(&untested), founders.to_string()];
    match &h.exact {
        Some(exact) => {
            let t = hap::total(exact);
            out.extend([t.ac.to_string(), num(t.af()), num(t.lower95()), num(t.upper95())]);
            out.extend(counts(exact));
            out.push(if t.ac > 0 {
                POPS[freq::strongest(exact)].to_string()
            } else {
                ".".to_string()
            });
            out.push(num(freq::heterogeneity_p(exact)));
        }
        None => out.extend(dots(11)),
    }
    if total > 0 && untested.len() == total {
        out.extend(dots(7));
    } else {
        let carrying = hap::total(&h.carrying);
        out.extend([carrying.ac.to_string(), num(carrying.af())]);
        out.extend(counts(&h.carrying));
    }
    out
}

/// `alignments.tsv`: where each allele was placed.
pub fn alignments_table(build: &str) -> Table<Placed> {
    let build = build.to_string();
    let placement = |cell: fn(&Placement) -> String| move |p: &Placed| p.placement.as_ref().map_or(".".to_string(), cell);
    Table::new()
        .col("allele", |p: &Placed| p.allele.clone())
        .col("gene", |p: &Placed| p.gene.clone())
        .col("build", move |_| build.clone())
        .col("contig", |p: &Placed| p.contig().to_string())
        .col("ref_start", placement(|a| a.span.0.to_string()))
        .col("ref_end", placement(|a| a.span.1.to_string()))
        .col("strand", placement(|a| if a.reverse { "-" } else { "+" }.to_string()))
        .col("identity", placement(|a| format!("{:.4}", a.identity)))
        .col("score", placement(|a| a.score.to_string()))
        .col(
            "best_other_gene",
            placement(|a| a.rival.as_ref().map_or(".".to_string(), |(gene, _)| gene.clone())),
        )
        .col(
            "best_other_score",
            placement(|a| a.rival.as_ref().map_or(".".to_string(), |(_, s)| s.to_string())),
        )
        .col("shares_locus_with", |p: &Placed| p.shared())
        .col("status", |p: &Placed| p.status.as_str().to_string())
        .col("missing_5prime", placement(|a| a.missing.0.to_string()))
        .col("missing_3prime", placement(|a| a.missing.1.to_string()))
}

/// `alignments.txt`: stacked alignments, one block per gene, to read by eye.
pub fn alignments_text(build: &str, windows: &[Window], placed: &[Placed]) -> Result<String> {
    let mut out = String::new();
    writeln!(out, "# igpop alignments, one block per gene, shown 5'->3' in gene orientation.")?;
    writeln!(
        out,
        "# '.' = same as reference, letter = mismatch, '-' = reference base missing in the allele,"
    )?;
    writeln!(
        out,
        "# blank = not covered by the allele. Inserted bases are listed under each block as"
    )?;
    writeln!(out, "# <reference position before the insertion>+<bases on the reference strand>.")?;
    writeln!(
        out,
        "# Coordinates above each chunk are the reference positions of its first and last column.\n"
    )?;
    let mut order: Vec<usize> = (0..windows.len()).collect();
    order.sort_by(|x, y| windows[*x].coord.gene.cmp(&windows[*y].coord.gene));
    for index in order {
        let members: Vec<(&Placed, &Placement)> = placed
            .iter()
            .filter_map(|p| p.placement.as_ref().filter(|a| a.window == index).map(|a| (p, a)))
            .collect();
        let Some((_, first)) = members.first() else { continue };
        let w = &windows[index];
        let reverse = first.reverse;
        let title = format!(
            "{}  {build} {}:{}-{}  strand {}  {} allele(s)",
            w.coord.gene,
            w.coord.contig,
            w.coord.start,
            w.coord.end,
            if reverse { '-' } else { '+' },
            members.len(),
        );
        let entries: Vec<render::Entry> = members
            .iter()
            .map(|(p, a)| render::Entry {
                name: &p.allele,
                row: &a.row,
                note: format!(
                    "identity {:.4}  status {}  variants {}  insertions {}{}",
                    a.identity,
                    p.status.as_str(),
                    p.variants.len(),
                    if a.insertions.is_empty() {
                        "none".to_string()
                    } else {
                        a.insertions.join(",")
                    },
                    if a.reverse == reverse {
                        ""
                    } else {
                        "  (aligned on the other strand)"
                    },
                ),
            })
            .collect();
        write!(out, "{}", render::gene_block(&title, &w.seq, w.start, reverse, &entries, 100))?;
    }
    Ok(out)
}

/// `score_calibration.tsv`: the counts of this run and the scores in use.
pub fn calibration_table(model: &ScoreModel) -> Result<String> {
    let mut out = String::new();
    writeln!(
        out,
        "support\tarbitrary_changes\tarbitrary_fraction\ttrusted_variants\ttrusted_fraction\tscore\tpooled_group\tscore_source"
    )?;
    let (na, nt): (u64, u64) = (model.arbitrary.iter().sum(), model.trusted.iter().sum());
    let fraction = |x: u64, n: u64| {
        if n == 0 {
            ".".to_string()
        } else {
            format!("{:.5}", x as f64 / n as f64)
        }
    };
    for i in 0..score::BINS {
        let (value, group) = match &model.calibration {
            Some(c) => (format!("{:.3}", c.score[i]), c.group[i].to_string()),
            None => (".".to_string(), ".".to_string()),
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{value}\t{group}\t{}",
            score::LABELS[i],
            model.arbitrary[i],
            fraction(model.arbitrary[i], na),
            model.trusted[i],
            fraction(model.trusted[i], nt),
            model.source.as_str()
        )?;
    }
    Ok(out)
}

/// `dbsnp_datasets.tsv`: the datasets relayed by dbSNP and their use.
pub fn datasets_table(ev: &Evidence) -> Result<String> {
    let mut out = String::new();
    writeln!(out, "dataset\trecords\tsmallest_nonzero_AF\timplied_chromosomes\tused\treason")?;
    for d in &ev.dbsnp.datasets {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}",
            d.name,
            d.records,
            num(d.smallest_af),
            d.chromosomes,
            if d.used { "yes" } else { "no" },
            d.reason
        )?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_and_rows_have_the_same_columns() {
        let table: Table<u32> = Table::new()
            .col("n", |n: &u32| n.to_string())
            .block(names(&["double", "square"]), |n: &u32| {
                vec![(2 * n).to_string(), (n * n).to_string()]
            });
        assert_eq!(table.header(), "n\tdouble\tsquare");
        assert_eq!(table.row(&3), "3\t6\t9");
    }

    #[test]
    #[should_panic(expected = "must fill all its columns")]
    fn a_block_that_fills_too_few_columns_is_caught() {
        let table: Table<u32> = Table::new().block(names(&["a", "b"]), |_| vec![String::new()]);
        table.row(&1);
    }

    #[test]
    fn numbers() {
        assert_eq!(num(f64::NAN), ".");
        assert_eq!(num(0.0), "0.000000");
        assert_eq!(num(0.00001234), "1.234e-5");
    }
}
