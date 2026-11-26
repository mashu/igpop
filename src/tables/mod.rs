//! Step 4: the output tables. Every table is declared as a list of columns,
//! each a name with the function that fills it, so a header and its rows
//! cannot drift apart. The meaning of every column is in docs/src/outputs.md.
//!
//! - `alleles`: one row per allele (`alleles.tsv`, `final.tsv`);
//! - `variants`: one row per difference of an allele (`variants.tsv`);
//! - `populations`: what each population says about each difference
//!   (`populations.tsv`), the quality of each position in each dataset
//!   (`site_quality.tsv`), and who the populations are (`datasets.tsv`);
//! - `alignments`: where each allele was placed, as a table and as text;
//! - `calibration`: the score per level of support and the haplotype
//!   check's calibration.

mod alignments;
mod alleles;
mod calibration;
mod populations;
mod variants;

pub use alignments::{alignments_legend, alignments_table, alignments_text, PlacedRow};
pub use alleles::alleles_table;
pub use calibration::{haplotype_table, levels_table};
pub use populations::{datasets_table, populations_table, quality_table, DatasetRow, PopulationRow, QualityRow, Site};
pub use variants::{variants_table, VariantRow};

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
    pub(crate) fn col(mut self, name: &str, cell: impl Fn(&T) -> String + 'static) -> Self {
        self.names.push(name.to_string());
        self.blocks.push((1, Box::new(move |row| vec![cell(row)])));
        self
    }

    /// Several columns filled together (one source, one case distinction).
    pub(crate) fn block(mut self, names: Vec<String>, cells: impl Fn(&T) -> Vec<String> + 'static) -> Self {
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

pub(crate) fn dots(n: usize) -> Vec<String> {
    vec![".".to_string(); n]
}

pub(crate) fn list(ids: &[String]) -> String {
    if ids.is_empty() {
        ".".to_string()
    } else {
        ids.join(";")
    }
}

pub(crate) fn names(columns: &[&str]) -> Vec<String> {
    columns.iter().map(|c| c.to_string()).collect()
}

/// A cell that may have no value.
pub(crate) fn or_dot<T>(value: Option<T>, text: impl Fn(T) -> String) -> String {
    value.map_or(".".to_string(), text)
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
