import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


class ImageStackViewer(tk.Tk):
    def __init__(self, image_stack):
        super().__init__()
        self.title("Image Stack Viewer")

        self.image_stack = image_stack
        self.num_images, self.image_size_x, self.image_size_y = image_stack.shape

        self.current_image_index = 0

        self.create_widgets()

    def create_widgets(self):
        # Matplotlib figure and canvas
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.canvas.draw()

        # Slider to navigate through images
        self.slider_label = ttk.Label(self, text="Image Index:")
        self.slider_label.pack(pady=10)
        self.image_slider = ttk.Scale(self, from_=0, to=self.num_images - 1, orient=tk.HORIZONTAL,
                                      command=self.update_image)
        self.image_slider.set(0)
        self.image_slider.pack(fill=tk.X, padx=10)

        # Navigation buttons
        self.prev_button = ttk.Button(self, text="Previous", command=self.show_previous_image)
        self.prev_button.pack(side=tk.LEFT, padx=10)
        self.next_button = ttk.Button(self, text="Next", command=self.show_next_image)
        self.next_button.pack(side=tk.RIGHT, padx=10)

        # Initialize the first image
        self.show_image()

    def show_image(self):
        self.ax.clear()
        self.ax.imshow(self.image_stack[self.current_image_index], cmap='gray')
        self.ax.set_title(f"Image {self.current_image_index + 1}/{self.num_images}")
        self.canvas.draw()

    def update_image(self, event):
        self.current_image_index = int(self.image_slider.get())
        self.show_image()

    def show_previous_image(self):
        self.current_image_index = (self.current_image_index - 1) % self.num_images
        self.image_slider.set(self.current_image_index)
        self.show_image()

    def show_next_image(self):
        self.current_image_index = (self.current_image_index + 1) % self.num_images
        self.image_slider.set(self.current_image_index)
        self.show_image()


