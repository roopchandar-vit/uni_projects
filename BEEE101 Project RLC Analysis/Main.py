import tkinter as tk
from tkinter import ttk
from math import pi, atan2, degrees
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image, ImageTk

# IMPEDENCE CALCULATIONS
def j(x):
    return 1j * x

def Z_series(R, L, C, omega):
    if C == 0 or omega == 0:
        return R + j(omega * L)
    return R + j(omega*L - 1/(omega*C))

def Z_parallel(R, L, C, omega):
    Y = 0
    if R != 0:
        Y += 1.0 / R
    if L != 0:
        Y += 1.0 / j(omega*L)
    if C != 0:
        Y += j(omega*C)
    if Y == 0:
        return complex(np.inf)
    return 1.0 / Y

def Z_R_series_parallel_LC(R, L, C, omega):
    if L == 0 and C == 0:
        return R
    ZL = j(omega*L) if L != 0 else complex(np.inf)
    ZC = 1.0 / (j(omega*C)) if (C != 0 and omega != 0) else complex(np.inf)
    if ZL == complex(np.inf):
        Zpar = ZC
    elif ZC == complex(np.inf):
        Zpar = ZL
    else:
        Zpar = 1.0 / (1.0/ZL + 1.0/ZC)
    return R + Zpar

def mag_angle(z):
    return abs(z), degrees(atan2(z.imag, z.real))

# MAIN APPLICATION
class PhasorApp:
    def __init__(self, root):
        self.root = root
        root.title("Interactive AC Circuit")
        root.geometry("1100x700")
        self.style = ttk.Style(root)
        try:
            self.style.theme_use("clam")
        except:
            pass
        img1 = "RLC-series.png"
        img2 = "parallel-RLC.png"
        img3 = "parallel-LC-series-R.png"
        o_img1=Image.open(img1).resize((300,300))
        o_img2=Image.open(img2).resize((300,300))
        o_img3=Image.open(img3).resize((300,300))
        self.tk_img1 = ImageTk.PhotoImage(o_img1)
        self.tk_img2 = ImageTk.PhotoImage(o_img2)
        self.tk_img3 = ImageTk.PhotoImage(o_img3)

        top = ttk.Frame(root, padding=8)
        top.pack(side="top", fill="x")
        ttk.Label(top, text="AC Circuit Simulator", font=("Segoe UI", 48, "bold")).pack()

        body = ttk.Frame(root, padding=8)
        body.pack(fill="both", expand=True)

        # CIRCUIT CHOICE
        left = ttk.Frame(body, padding=8)
        left.pack(fill="y")
        ttk.Label(left, text="Select Circuit", font=("Segoe UI", 24, "bold")).pack(pady=(0,8))

        ttk.Button(left, image=self.tk_img1, text="LCR Series", compound="top", width=25,
                   command=lambda: self.open_circuit("series")).pack(side="left", pady=10)
        ttk.Button(left, image=self.tk_img2, text="LCR Parallel", compound="top", width=25,
                   command=lambda: self.open_circuit("parallel")).pack(side="left", pady=10)
        ttk.Button(left, image=self.tk_img3, text="R in series with (L‖C)", compound="top", width=25,
                   command=lambda: self.open_circuit("r_par_lc")).pack(side="left", pady=10)

        ttk.Separator(left, orient="horizontal").pack(fill="x", pady=8)

        self.right = ttk.Frame(body, padding=8, relief="flat")
        self.right.pack(side="left", fill="both", expand=True)
        self.circuit_windows = {}

    def open_circuit(self, circuit_type):
        if circuit_type in self.circuit_windows and self.circuit_windows[circuit_type].winfo_exists():
            self.circuit_windows[circuit_type].lift()
            return

        win = tk.Toplevel(self.root)
        win.title(f"Circuit: {circuit_type.upper()}")
        win.geometry("1000x700")
        self.circuit_windows[circuit_type] = win

        ctrl = ttk.Frame(win, padding=6)
        ctrl.pack(side="top", fill="x")

        R_var = tk.DoubleVar(value=10.0)
        L_var = tk.DoubleVar(value=0.01)
        C_var = tk.DoubleVar(value=1e-6)
        V_var = tk.DoubleVar(value=230.0)
        f_var = tk.DoubleVar(value=50.0)
        show_component = tk.BooleanVar(value=True)

        def make_scale(label, var, frm, row, rng):
            ttk.Label(frm, text=label).grid(row=row, column=0, sticky="w", padx=4)
            scale = ttk.Scale(frm, variable=var, from_=rng[0], to=rng[1], orient="horizontal")
            scale.grid(row=row, column=1, sticky="ew", padx=4)
            ttk.Entry(frm, textvariable=var, width=10).grid(row=row, column=2, padx=6)
            return scale

        param_frame = ttk.Frame(ctrl)
        param_frame.pack(side="left", fill="x", expand=True, padx=8)
        param_frame.columnconfigure(1, weight=1)

        make_scale("R (Ω):", R_var, param_frame, 0, (0.1, 1000))
        make_scale("L (H):", L_var, param_frame, 1, (1e-6, 0.5))
        make_scale("C (F):", C_var, param_frame, 2, (1e-9, 1e-3))
        make_scale("V (RMS):", V_var, param_frame, 3, (1, 1000))
        make_scale("Frequency (Hz):", f_var, param_frame, 4, (0.1, 1000))

        right_top = ttk.Frame(ctrl)
        right_top.pack(side="right", padx=8)
        ttk.Checkbutton(right_top, text="Show component phasors", variable=show_component).pack(pady=4)
        ttk.Button(right_top, text="Reset Defaults",
                   command=lambda: (R_var.set(10), L_var.set(0.01), C_var.set(1e-6), V_var.set(230), f_var.set(50))).pack(pady=4)

        info = ttk.Frame(win, padding=6)
        info.pack(side="top", fill="x")
        z_label = ttk.Label(info, text="Impedance: --  |  Current: --  |  Power (P/Q/S): -- / -- / --")
        z_label.pack(anchor="w", padx=6, pady=4)

        plots = ttk.Frame(win)
        plots.pack(fill="both", expand=True)
        fig = Figure(figsize=(8,6), dpi=100)
        ax_time = fig.add_subplot(211)
        ax_phasor = fig.add_subplot(212)
        fig.tight_layout(pad=3.0)

        v_line, = ax_time.plot([], [], label="v(t)", color='tab:orange')
        i_line, = ax_time.plot([], [], label="i(t)", color='tab:blue')
        p_line, = ax_time.plot([], [], label="p(t)", color='tab:green')
        ax_time.legend()
        ax_time.set_xlabel("Time (s)")
        ax_time.set_ylabel("Amplitude")

        canvas = FigureCanvasTkAgg(fig, master=plots)
        canvas.get_tk_widget().pack(fill="both", expand=True)

        running = {"on": True}
        fps = 15

        def compute_and_update():
            R, L, C = R_var.get(), L_var.get(), C_var.get()
            V_rms, f = V_var.get(), f_var.get()
            omega = 2*pi*f
            V_peak = V_rms*np.sqrt(2)

            if circuit_type == "series":
                Z = Z_series(R, L, C, omega)
            elif circuit_type == "parallel":
                Z = Z_parallel(R, L, C, omega)
            else:
                Z = Z_R_series_parallel_LC(R, L, C, omega)

            I = V_rms / Z
            Zm, Za = mag_angle(Z)
            Im, Ia = mag_angle(I)

            S = V_rms * I.conjugate()
            P, Q = S.real, S.imag
            z_label.config(text=f"Z={Zm:.3f}Ω ∠{Za:.2f}° | I={Im:.3f}A ∠{Ia:.2f}° | Real Power (P)/Reactive Power (Q)/Apparent Power(S)={P:.2f}/{Q:.2f}/{abs(S):.2f}")

            # time-domain update
            T = 1/max(f, 0.001)
            t = np.linspace(0, T, 1000)
            v_t = V_peak*np.sin(omega*t)
            i_t = Im*np.sqrt(2)*np.sin(omega*t + np.radians(Ia))
            p_t = v_t * i_t

            ax_time.clear()
            ax_time.plot(t, v_t, label="v(t)", color='tab:orange')
            ax_time.plot(t, i_t, label="i(t)", color='tab:blue')
            ax_time.plot(t, p_t, label="p(t)", color='tab:green')
            ax_time.legend(loc="upper right")
            ax_time.set_xlabel("Time (s)")
            ax_time.set_ylabel("Amplitude")

            # phasor plot
            ax_phasor.clear()
            ax_phasor.grid(True)
            ax_phasor.set_aspect('equal')
            Ix, Iy = np.cos(0), np.sin(0)
            Vx, Vy = np.cos(np.radians(Za)), np.sin(np.radians(Za))
            ax_phasor.arrow(0, 0, Vx, Vy, head_width=0.05, color='tab:orange', lw=3, length_includes_head=True)
            ax_phasor.arrow(0, 0, Ix, Iy, head_width=0.05, color='tab:blue', lw=3, length_includes_head=True)
            ax_phasor.text(Vx*1.1, Vy*1.1, "V", color='tab:orange', fontsize=12)
            ax_phasor.text(Ix*1.1, Iy*1.1, "I", color='tab:blue', fontsize=12)
            ax_phasor.set_xlim(-1.5, 1.5)
            ax_phasor.set_ylim(-1.5, 1.5)
            ax_phasor.set_title("Phasor Diagram (I reference)")
            canvas.draw_idle()

        def loop():
            if not running["on"]:
                return
            try:
                compute_and_update()
            except Exception as e:
                print("Update error:", e)
            win.after(int(1000/fps), loop)

        loop()

        def on_close():
            running["on"] = False
            win.destroy()
        win.protocol("WM_DELETE_WINDOW", on_close)

root = tk.Tk()
app = PhasorApp(root)
root.mainloop()
