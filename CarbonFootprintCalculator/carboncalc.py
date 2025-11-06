import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
from matplotlib import pyplot as plt

# Emission factors (kg CO₂ per unit)
FACTORS = {"travel": 0.21, "electricity": 0.475, "diet_nonveg": 2.5, "diet_veg": 1.0}

# Main calculation function
def calculate():
    try:
        travel = float(travel_entry.get())
        electricity = float(electricity_entry.get())
        diet = diet_choice.get()

        travel_em = travel * FACTORS["travel"] * 365
        elec_em = electricity * FACTORS["electricity"] * 365
        diet_factor = FACTORS["diet_nonveg"] if diet == "Non-Vegetarian" else FACTORS["diet_veg"]
        diet_em = diet_factor * 3 * 365 / 7

        total = travel_em + elec_em + diet_em
        tons = total / 1000

        # Display result
        result_label.config(text=f"Your annual carbon footprint: {tons:.2f} tons CO₂/year")
        details_label.config(text=f"({round(total)} kg CO₂ per year)")

        # Category
        if tons < 3:
            category = "Low Footprint"
            color = "green"
        elif tons < 6:
            category = "Moderate Footprint"
            color = "orange"
        else:
            category = "High Footprint"
            color = "red"

        category_label.config(text=category, fg=color)

        # Suggestions
        suggestions.delete(0, tk.END)
        if travel > 20:
            suggestions.insert(tk.END, "Use public transport or carpool for long commutes.")
        elif travel > 5:
            suggestions.insert(tk.END, "Walk or cycle for short trips.")
        else:
            suggestions.insert(tk.END, "Great travel habits!")

        if electricity > 10:
            suggestions.insert(tk.END, "Switch to LED bulbs and efficient appliances.")
            suggestions.insert(tk.END, "Consider solar panels or green energy.")
        elif electricity > 5:
            suggestions.insert(tk.END, "Unplug chargers when not in use.")
        else:
            suggestions.insert(tk.END, "Electricity use is efficient!")

        if diet == "Non-Vegetarian":
            suggestions.insert(tk.END, "Try reducing meat consumption to a few times per week.")
        else:
            suggestions.insert(tk.END, "Vegetarian diet helps lower CO₂ impact.")

        # Pie Chart
        plt.figure(figsize=(4, 4))
        plt.pie(
            [travel_em, elec_em, diet_em],
            labels=["Travel", "Electricity", "Diet"],
            autopct="%1.1f%%",
            colors=["#16a34a", "#60a5fa", "#f59e0b"]
        )
        plt.title("CO₂ Emission Breakdown")
        plt.show()

    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers for all fields.")

def reset():
    travel_entry.delete(0, tk.END)
    travel_entry.insert(0, "10")
    electricity_entry.delete(0, tk.END)
    electricity_entry.insert(0, "6")
    diet_choice.set("Vegetarian")
    result_label.config(text="")
    category_label.config(text="")
    details_label.config(text="")
    suggestions.delete(0, tk.END)

# GUI Setup
root = tk.Tk()
root.title(" Carbon Footprint Calculator")
root.geometry("500x600")
root.config(bg="#0f1724")

# Title
title = tk.Label(root, text="Carbon Footprint Calculator", font=("Inter", 16, "bold"), fg="white", bg="#0f1724")
title.pack(pady=10)

# Input fields
frame = tk.Frame(root, bg="#0f1724")
frame.pack(pady=10)

tk.Label(frame, text="Daily Travel (km):", bg="#0f1724", fg="white").grid(row=0, column=0, sticky="w", pady=5)
travel_entry = tk.Entry(frame, width=20)
travel_entry.insert(0, "10")
travel_entry.grid(row=0, column=1, pady=5)

tk.Label(frame, text="Daily Electricity Use (kWh):", bg="#0f1724", fg="white").grid(row=1, column=0, sticky="w", pady=5)
electricity_entry = tk.Entry(frame, width=20)
electricity_entry.insert(0, "6")
electricity_entry.grid(row=1, column=1, pady=5)

tk.Label(frame, text="Diet Type:", bg="#0f1724", fg="white").grid(row=2, column=0, sticky="w", pady=5)
diet_choice = ttk.Combobox(frame, values=["Vegetarian", "Non-Vegetarian"], width=18)
diet_choice.set("Vegetarian")
diet_choice.grid(row=2, column=1, pady=5)

# Buttons
btn_frame = tk.Frame(root, bg="#0f1724")
btn_frame.pack(pady=10)

calc_btn = tk.Button(btn_frame, text="Calculate", command=calculate, bg="#16a34a", fg="white", width=12)
calc_btn.grid(row=0, column=0, padx=5)

reset_btn = tk.Button(btn_frame, text="Reset", command=reset, bg="#1e293b", fg="white", width=12)
reset_btn.grid(row=0, column=1, padx=5)

# Results
result_label = tk.Label(root, text="", font=("Inter", 12), bg="#0f1724", fg="white")
result_label.pack(pady=5)

details_label = tk.Label(root, text="", font=("Inter", 10), bg="#0f1724", fg="#9ca3af")
details_label.pack()

category_label = tk.Label(root, text="", font=("Inter", 12, "bold"), bg="#0f1724")
category_label.pack(pady=5)

# Suggestions list
tk.Label(root, text="Suggestions:", bg="#0f1724", fg="white", font=("Inter", 12, "bold")).pack(pady=5)
suggestions = tk.Listbox(root, width=60, height=6, bg="#1e293b", fg="white", border=0)
suggestions.pack(pady=5)

root.mainloop()