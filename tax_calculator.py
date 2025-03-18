import streamlit as st
import pandas as pd

def calculate_tax(ord_income, divcap_income, ss_income, adjustments, qbi, filing_status, year):
    """Calculates US federal income tax based on provided inputs."""

    # Tax Brackets (2024 & 2025) - Single
    single_2024 = {
        10: (0, 11600),
        12: (11601, 47150),
        22: (47151, 100525),
        24: (100526, 191950),
        32: (191951, 243725),
        35: (243726, 609350),
        37: (609351, float('inf'))
    }

    single_2025 = {
        10: (0, 11600),
        12: (11601, 47150),
        22: (47151, 100525),
        24: (100526, 191950),
        32: (191951, 243725),
        35: (243726, 609350),
        37: (609351, float('inf'))
    }

    # Tax Brackets (2024 & 2025) - Married Filing Jointly
    married_2024 = {
        10: (0, 23200),
        12: (23201, 94300),
        22: (94301, 191950),
        24: (191951, 383900),
        32: (383901, 487450),
        35: (487451, 731200),
        37: (731201, float('inf'))
    }

    married_2025 = {
        10: (0, 23200),
        12: (23201, 94300),
        22: (94301, 191950),
        24: (191951, 383900),
        32: (383901, 487450),
        35: (487451, 731200),
        37: (731201, float('inf'))
    }

    # Tax Brackets (2024 & 2025) - Head of Household
    head_2024 = {
        10: (0, 17400),
        12: (17401, 59850),
        22: (59851, 123750),
        24: (123751, 231250),
        32: (231251, 277550),
        35: (277551, 578125),
        37: (578126, float('inf'))
    }

    head_2025 = {
        10: (0, 17400),
        12: (17401, 59850),
        22: (59851, 123750),
        24: (123751, 231250),
        32: (231251, 277550),
        35: (277551, 578125),
        37: (578126, float('inf'))
    }

    if year == 2024:
        brackets = {
            "Single": single_2024,
            "Married": married_2024,
            "Head of Household": head_2024
        }
    else:
        brackets = {
            "Single": single_2025,
            "Married": married_2025,
            "Head of Household": head_2025
        }

    bracket = brackets[filing_status]

    # Standard Deduction (2024 & 2025)
    std_ded_2024 = {
        "Single": 14600,
        "Married": 29200,
        "Head of Household": 21900
    }

    std_ded_2025 = {
        "Single": 14600,
        "Married": 29200,
        "Head of Household": 21900
    }

    if year == 2024:
        std_ded = std_ded_2024[filing_status]
    else:
        std_ded = std_ded_2025[filing_status]

    # Calculations
    tot_ordinc_with_ss = ord_income + ss_income
    divcap = divcap_income
    ss_income = ss_income
    adjustments = adjustments
    qbi = qbi

    ordinc = tot_ordinc_with_ss - adjustments
    agi = ordinc - qbi
    ordinc_taxable = max(0, ordinc)
    divcap_taxable = max(0, divcap)
    total_taxable = ordinc_taxable + divcap_taxable
    total_ded = std_ded + qbi

    # Calculate Taxable Social Security
    if ss_income > 0:
        if filing_status == "Single":
            threshold = 25000
        elif filing_status == "Married":
            threshold = 32000
        else:  # Head of Household
            threshold = 25000

        ss_taxable = max(0, ss_income - threshold)
    else:
        ss_taxable = 0

    # Calculate Tax
    ord_tax = 0
    for rate, (lower, upper) in bracket.items():
        if ordinc_taxable > lower:
            taxable_in_bracket = min(ordinc_taxable, upper) - lower
            ord_tax += taxable_in_bracket * (rate / 100)

    divcap_tax = 0
    # Assuming qualified dividends and long-term capital gains are taxed at 0%, 15%, or 20%
    if divcap_taxable > 0:
        if divcap_taxable <= 47150:
            divcap_tax = divcap_taxable * 0.15
        else:
            divcap_tax = divcap_taxable * 0.20

    total_tax = ord_tax + divcap_tax
    avg_tax_rate = (total_tax / total_taxable) * 100 if total_taxable > 0 else 0
    cap_ded = 0  # Placeholder for capital loss deduction

    # Create Output Table
    data = {
        "Item": [
            "Ordinary Income",
            "Dividends and Capital Gains Income",
            "Social Security Income",
            "Taxable Social Security",
            "Total Income",
            "Income Adjustments",
            "Adj Ordinary Income",
            "AGI",
            "Standard Deduction",
            "QBI Deduction",
            "Total Deductions",
            "Taxable Ordinary Income",
            "Taxable Divcap",
            "Total Taxable Income",
            "Ordinary Income Tax",
            "Divcap Tax",
            "Total Tax",
            "Avg Tax Rate",
            "DivCap Deduction"
        ],
        "Value": [
            f"${tot_ordinc_with_ss:,.0f}",
            f"${divcap:,.0f}",
            f"${ss_income:,.0f}",
            f"${ss_taxable:,.0f}",
            f"${tot_ordinc_with_ss + divcap:,.0f}",
            f"${adjustments:,.0f}",
            f"${ordinc:,.0f}",
            f"${agi:,.0f}",
            f"${std_ded:,.0f}",
            f"${qbi:,.0f}",
            f"${total_ded:,.0f}",
            f"${ordinc_taxable:,.0f}",
            f"${divcap_taxable:,.0f}",
            f"${total_taxable:,.0f}",
            f"${ord_tax:,.0f}",
            f"${divcap_tax:,.0f}",
            f"${total_tax:,.0f}",
            f"{avg_tax_rate:.2f}%",
            f"{cap_ded:.2f}%"
        ],
        "Form Line": [
            "",
            "",
            "5a",
            "5b",
            "9",
            "10",
            "",
            "11",
            "12",
            "13",
            "14",
            "",
            "",
            "15",
            "",
            "",
            "16",
            "",
            ""
        ]
    }

    df = pd.DataFrame(data)
    df = df.set_index("Item")
    return df

# Streamlit App
st.title("US Federal Income Tax Calculator (2024/2025)")

# Inputs
ord_income = st.number_input("Ordinary Income", value=75000)
divcap_income = st.number_input("Dividends and Capital Gains Income", value=10000)
ss_income = st.number_input("Social Security Income", value=0)
adjustments = st.number_input("Income Adjustments (e.g., HSA)", value=0)
qbi = st.number_input("Qualified Business Income (QBI) Deduction", value=0)

filing_status = st.selectbox("Filing Status", ["Single", "Married", "Head of Household"])
year = st.selectbox("Tax Year", [2024, 2025])

# Calculate Tax
if st.button("Calculate"):
    tax_df = calculate_tax(ord_income, divcap_income, ss_income, adjustments, qbi, filing_status, year)
    st.dataframe(tax_df)
    