import argparse
import re
from pathlib import Path

from astropy.table import Table


TABLE_DIR = Path(__file__).resolve().parent.parent / "tables"

table_primary = Table.read(TABLE_DIR / "table_primary.csv")
table_extension = Table.read(TABLE_DIR / "table_extension.csv")
table_footnotes = Table.read(TABLE_DIR / "table_footnotes.csv")
table_legend = Table.read(TABLE_DIR / "table_legend.csv")

# These CSV files currently contain a UTF-8 byte-order mark before "key".
table_footnotes.rename_column(table_footnotes.colnames[0], "key")
table_legend.rename_column(table_legend.colnames[0], "key")


# Codes present in the keyword tables but not currently defined in
# table_legend.csv.
LEGACY_LEGEND = {
    "M": "Keyword mandatory.",
    "RC": "Keyword recommended.",
    "SU": "Contact Phase 3 support to determine whether the keyword is required.",
}


def _get_legend_lookup(legend_table):
    legend = {
        str(key).strip(): str(meaning).strip()
        for key, meaning in zip(legend_table["key"], legend_table["meaning"])
    }
    return {**LEGACY_LEGEND, **legend}


def _describe_requirement(requirement, legend_lookup):
    """Expand a code such as ``MA65`` into a meaning and footnote number."""
    descriptions = []
    for code in requirement.split():
        match = re.fullmatch(r"([A-Z]+)(\d*)", code)
        if match is None:
            descriptions.append({"code": code, "meaning": "Unknown requirement code"})
            continue

        legend_key, number = match.groups()
        description = {
            "code": code,
            "legend_key": legend_key,
            "meaning": legend_lookup.get(
                legend_key, f"Unknown requirement code {legend_key!r}"
            ),
        }
        if number:
            description["footnote"] = int(number)
        descriptions.append(description)
    return descriptions


def _format_requirement_detail(detail):
    """Format a decoded requirement for the human-readable report."""
    meaning = detail["meaning"].rstrip(".")
    if "footnote" in detail:
        # The footnote marker printed below makes this clause redundant.
        meaning = meaning.split(", refer to the table footnotes", 1)[0]
        meaning += f" [footnote {detail['footnote']}]"
    return meaning


def _get_keyword_entries(product_name, table, legend_lookup):
    """Return the applicable keyword entries for one product column."""
    if product_name not in table.colnames:
        available = ", ".join(table.colnames[2:])
        raise ValueError(
            f"Unknown product type {product_name!r}. Available types: {available}"
        )

    entries = []
    for name, requirement in zip(table["Column Name"], table[product_name]):
        requirement = str(requirement).strip()
        if requirement != "NA":
            entries.append(
                {
                    "keyword": str(name).strip(),
                    "requirement": requirement,
                    "requirement_details": _describe_requirement(
                        requirement, legend_lookup
                    ),
                }
            )
    return entries


def _get_footnotes(entries, footnote_table):
    """Return footnotes referenced by requirement codes such as MA65 or M33."""
    referenced_numbers = {
        int(number)
        for entry in entries
        for number in re.findall(r"\d+", entry["requirement"])
    }
    footnote_lookup = {
        int(key): str(text).strip()
        for key, text in zip(footnote_table["key"], footnote_table["footnote"])
    }

    return [
        {"number": number, "text": footnote_lookup[number]}
        for number in sorted(referenced_numbers)
        if number in footnote_lookup
    ]


def get_values(
    product_name,
    primary_table=table_primary,
    extension_table=table_extension,
    footnote_table=table_footnotes,
    legend_table=table_legend,
    *,
    print_summary=True,
):
    """Return keyword requirements and referenced footnotes for a product type.

    Entries whose table value is ``NA`` are omitted. Requirement codes are kept
    in the result so that mandatory, optional, recommended, and other entries
    remain distinguishable.
    """
    legend_lookup = _get_legend_lookup(legend_table)
    primary = _get_keyword_entries(product_name, primary_table, legend_lookup)
    extension = _get_keyword_entries(product_name, extension_table, legend_lookup)
    footnotes = _get_footnotes(primary + extension, footnote_table)

    result = {
        "product_name": product_name,
        "primary": primary,
        "extension": extension,
        "footnotes": footnotes,
        "legend": legend_lookup,
    }

    if print_summary:
        for heading, entries in (
            ("Primary HDU", primary),
            ("Extension HDU", extension),
        ):
            print("-" * 30)
            print(f"{product_name} - {heading}")
            print("-" * 30)
            for entry in entries:
                descriptions = [
                    _format_requirement_detail(detail)
                    for detail in entry["requirement_details"]
                ]
                print(f"{entry['keyword']}: {'; '.join(descriptions)}")
            print()

        print("-" * 30)
        print(f"{product_name} - Referenced footnotes")
        print("-" * 30)
        for footnote in footnotes:
            print(f"{footnote['number']}: {footnote['text']}")
        print()

    return result


def main():
    product_types = table_primary.colnames[2:]
    parser = argparse.ArgumentParser(
        description="Show the header keyword requirements for an ESO product type."
    )
    parser.add_argument(
        "product_name",
        choices=product_types,
        help="Product type to inspect.",
    )
    args = parser.parse_args()
    get_values(args.product_name)


if __name__ == "__main__":
    main()
