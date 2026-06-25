import re
import time
import requests
from pathlib import Path
import pdfplumber
import pandas as pd
from datetime import datetime


def extract_compounds_from_pdf(pdf_path):
    pdf_name = Path(pdf_path).name
    print(f"Обработка: {pdf_name}")

    compounds = []
    found_names = set()

    with pdfplumber.open(pdf_path) as pdf:
        full_text = ""
        for page in pdf.pages:
            full_text += page.extract_text() + "\n"

        lines = full_text.split("\n")

        name_patterns = [
            r"\b([A-Z][a-z]?[A-Z]?[a-z]?[\d\-]*[A-Z][a-z]?[\d\-]*[A-Z][a-z]?[\d\-]*[A-Z][a-z]?[\d\-]*)\b",
            r"\b([A-Za-z]+[\-]?[A-Za-z]+[\-]?[A-Za-z]+[\-]?[A-Za-z]+(?:amide|amine|acid|one|ol|ene|ane|ine|yl|ate|ide))\b",
            r"\b([A-Z][a-z]+[A-Z][a-z]+[A-Z][a-z]+)\b",
            r"\b([Nn]\-[A-Za-z0-9\-\(\)]+[A-Za-z0-9\-\(\)]+[A-Za-z0-9\-\(\)]+)\b",
            r"\b(\d+[\-]\([A-Za-z0-9\-]+\)[\-][A-Za-z0-9\-\(\)\[\]]+)\b",
            r"\b([A-Z][a-z]?\d*[A-Z]?[a-z]?\d*[A-Z]?[a-z]?\d*[A-Z]?[a-z]?\d*)\b",
        ]

        activity_patterns = [
            (r"IC50\s*[:=]\s*([\d.]+)\s*[nµ]?M", "IC50"),
            (r"Ki\s*[:=]\s*([\d.]+)\s*nM", "Ki"),
            (r"EC50\s*[:=]\s*([\d.]+)\s*[nµ]?M", "EC50"),
            (r"pIC50\s*[:=]\s*([\d.]+)", "pIC50"),
            (r"pKi\s*[:=]\s*([\d.]+)", "pKi"),
        ]

        for line in lines:
            for pattern in name_patterns:
                matches = re.findall(pattern, line)
                for name in matches:
                    name = name.strip()
                    if len(name) > 5 and name not in found_names:
                        found_names.add(name)

                        activities = []
                        for act_pattern, label in activity_patterns:
                            act_match = re.search(act_pattern, line, re.IGNORECASE)
                            if act_match:
                                activities.append(f"{label}={act_match.group(1)}")

                        if not activities:
                            current_idx = lines.index(line)
                            for offset in range(1, 4):
                                if current_idx + offset < len(lines):
                                    next_line = lines[current_idx + offset]
                                    for act_pattern, label in activity_patterns:
                                        act_match = re.search(
                                            act_pattern, next_line, re.IGNORECASE
                                        )
                                        if act_match:
                                            activities.append(
                                                f"{label}={act_match.group(1)}"
                                            )
                                            break
                                    if activities:
                                        break

                        compounds.append(
                            {
                                "name": name,
                                "activities": (
                                    "|".join(activities) if activities else ""
                                ),
                                "source": pdf_name,
                                "context": (
                                    line[:150] + "..." if len(line) > 150 else line
                                ),
                            }
                        )

    return compounds


def remove_duplicates(compounds):
    """Удаляет повторяющиеся соединения по названию"""
    print("\n🔍 Удаление дубликатов...")

    seen = set()
    unique = []

    for comp in compounds:
        name_lower = comp["name"].lower()
        if name_lower not in seen:
            seen.add(name_lower)
            unique.append(comp)

    removed = len(compounds) - len(unique)
    if removed > 0:
        print(f"  ✅ Удалено {removed} дубликатов")
    else:
        print(f"  ✅ Дубликатов не найдено")

    return unique


def check_in_pubchem(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/MolecularFormula/JSON"
    response = requests.get(url, timeout=10)

    if response.status_code == 200:
        data = response.json()
        properties = data.get("PropertyTable", {}).get("Properties", [])
        if properties:
            return True

    if "-" in name:
        parts = name.split("-")
        for part in parts[-2:]:
            if len(part) > 3:
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{part}/property/MolecularFormula/JSON"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    properties = data.get("PropertyTable", {}).get("Properties", [])
                    if properties:
                        return True

    return False


def validate_compounds_in_pubchem(compounds):
    print("\n🔍 Проверка соединений в PubChem...")

    validated = []
    total = len(compounds)

    for i, comp in enumerate(compounds, 1):
        print(f"  [{i}/{total}] Проверка: {comp['name'][:50]}...")

        if check_in_pubchem(comp["name"]):
            print(f"  Найдено в PubChem")
            comp["in_pubchem"] = True
            validated.append(comp)
        else:
            print(f"  Не найдено в PubChem - удаляем")

        time.sleep(0.5)

    print(f"После проверки осталось {len(validated)} соединений из {total}")
    return validated


def parse_all_pdfs():
    all_compounds = []
    pdf_files = list(Path(".").glob("*.pdf"))

    print(f"Найдено {len(pdf_files)} PDF файлов")

    for pdf_file in pdf_files:
        compounds = extract_compounds_from_pdf(pdf_file)
        print(f" найдено {len(compounds)} соединений")
        for comp in compounds[:3]:
            print(f"    - {comp['name'][:60]}...")
            if comp["activities"]:
                print(f"      Активности: {comp['activities']}")
        all_compounds.extend(compounds)

    return all_compounds


def save_to_csv(compounds, output_file="compound_database.csv"):
    df = pd.DataFrame(compounds)

    df = df.rename(
        columns={
            "name": "Name",
            "activities": "Linked_BioAssays",
            "source": "Data_Source",
            "context": "Context",
            "in_pubchem": "In_PubChem",
        }
    )

    df["cid"] = [f"C{i+1}" for i in range(len(df))]
    df["Synonyms"] = df["Name"]
    df["Create_Date"] = datetime.now().strftime("%Y%m%d")

    columns = [
        "cid",
        "Name",
        "Synonyms",
        "In_PubChem",
        "Linked_BioAssays",
        "Data_Source",
        "Context",
        "Create_Date",
    ]

    for col in columns:
        if col not in df.columns:
            df[col] = ""

    df = df[columns]
    df.to_csv(output_file, index=False, encoding="utf-8-sig")

    print(f"База данных сохранена в {output_file}")
    print(f" Всего соединений: {len(df)}")
    print(f"  Найдено в PubChem: {df[df['In_PubChem'] == True].shape[0]}")
    print(f"  Не найдено в PubChem: {df[df['In_PubChem'] == False].shape[0]}")

    print("Примеры соединений с активностями:")
    sample = df[df["Linked_BioAssays"] != ""].head(10)
    if not sample.empty:
        print(
            sample[["cid", "Name", "Linked_BioAssays", "Data_Source"]]
            .head(10)
            .to_string(index=False)
        )
    else:
        print("  (Активности не найдены)")


def main():
    print("ИЗВЛЕЧЕНИЕ СОЕДИНЕНИЙ ИЗ PDF С ПРОВЕРКОЙ В PUBCHEM")

    compounds = parse_all_pdfs()

    # Удаляем дубликаты
    compounds = remove_duplicates(compounds)

    validated = validate_compounds_in_pubchem(compounds)
    save_to_csv(validated)


if __name__ == "__main__":
    main()
