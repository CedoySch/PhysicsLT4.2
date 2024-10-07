import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QLabel, QLineEdit, QGridLayout, QMessageBox, QFrame
from PyQt5.QtCore import QTimer
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class ElasticCollisionSimulation(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Настройка и симуляция упругого соударения')

        self.layout = QVBoxLayout()

        self.grid_layout = QGridLayout()

        # Настройки для тела 1
        self.mass1_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Масса тела 1 (кг):'), 0, 0)
        self.grid_layout.addWidget(self.mass1_input, 0, 1)

        self.v0_1_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Начальная скорость тела 1 (м/с):'), 1, 0)
        self.grid_layout.addWidget(self.v0_1_input, 1, 1)

        self.angle1_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Угол направления тела 1 (градусы):'), 2, 0)
        self.grid_layout.addWidget(self.angle1_input, 2, 1)

        self.radius1_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Радиус тела 1 (м):'), 3, 0)
        self.grid_layout.addWidget(self.radius1_input, 3, 1)

        # Линия-разделитель
        line1 = QFrame()
        line1.setFrameShape(QFrame.HLine)
        line1.setFrameShadow(QFrame.Sunken)
        self.grid_layout.addWidget(line1, 4, 0, 1, 2)

        # Настройки для тела 2
        self.mass2_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Масса тела 2 (кг):'), 5, 0)
        self.grid_layout.addWidget(self.mass2_input, 5, 1)

        self.v0_2_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Начальная скорость тела 2 (м/с):'), 6, 0)
        self.grid_layout.addWidget(self.v0_2_input, 6, 1)

        self.angle2_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Угол направления тела 2 (градусы):'), 7, 0)
        self.grid_layout.addWidget(self.angle2_input, 7, 1)

        self.radius2_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Радиус тела 2 (м):'), 8, 0)
        self.grid_layout.addWidget(self.radius2_input, 8, 1)

        # Линия-разделитель
        line2 = QFrame()
        line2.setFrameShape(QFrame.HLine)
        line2.setFrameShadow(QFrame.Sunken)
        self.grid_layout.addWidget(line2, 9, 0, 1, 2)

        # Настройка ширины и высоты контейнера
        self.width_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Ширина контейнера (м):'), 10, 0)
        self.grid_layout.addWidget(self.width_input, 10, 1)

        self.height_input = QLineEdit()
        self.grid_layout.addWidget(QLabel('Высота контейнера (м):'), 11, 0)
        self.grid_layout.addWidget(self.height_input, 11, 1)

        self.layout.addLayout(self.grid_layout)

        self.canvas = FigureCanvas(plt.figure())
        self.ax = self.canvas.figure.add_subplot(111)
        self.ax.set_xlim(0, 40)
        self.ax.set_ylim(0, 10)
        self.ax.set_aspect('equal')

        self.layout.addWidget(self.canvas)

        self.start_button = QPushButton('Запустить симуляцию')
        self.start_button.clicked.connect(self.start_simulation)
        self.layout.addWidget(self.start_button)

        self.setLayout(self.layout)

        self.timer = QTimer()
        self.timer.timeout.connect(self.update_simulation)

        self.body1 = None
        self.body2 = None
        self.initial_energy = 0  # Изначальная энергия

    def show_error(self, message):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.setWindowTitle('Ошибка')
        msg_box.setText(message)
        msg_box.exec_()

    def start_simulation(self):
        try:
            self.m1 = float(self.mass1_input.text())
            self.m2 = float(self.mass2_input.text())
            v0_1 = float(self.v0_1_input.text())
            angle1 = float(self.angle1_input.text())
            v0_2 = float(self.v0_2_input.text())
            angle2 = float(self.angle2_input.text())
            self.r1 = float(self.radius1_input.text())
            self.r2 = float(self.radius2_input.text())
            self.width = float(self.width_input.text())
            self.height = float(self.height_input.text())

            # Проверка ограничений
            if self.m1 <= 0 or self.m2 <= 0:
                raise ValueError("Масса должна быть положительным числом.")
            if self.r1 <= 0 or self.r2 <= 0:
                raise ValueError("Радиус должен быть положительным числом.")
            if v0_1 == 0 or v0_2 == 0:
                raise ValueError("Скорость не может быть равна нулю.")
            if self.width <= self.r1 + self.r2 or self.height <= self.r1 + self.r2:
                raise ValueError("Ширина и высота контейнера должны быть больше суммы радиусов тел.")
        except ValueError:
            self.show_error("Введите корректные числовые значения")
            return

        angle1 = np.radians(angle1)
        angle2 = np.radians(angle2)

        self.vx1 = v0_1 * np.cos(angle1)
        self.vy1 = v0_1 * np.sin(angle1)
        self.vx2 = v0_2 * np.cos(angle2)
        self.vy2 = v0_2 * np.sin(angle2)

        self.x1, self.y1 = 2, 5
        self.x2, self.y2 = 7, 5

        self.ax.clear()
        self.ax.set_xlim(0, self.width)
        self.ax.set_ylim(0, self.height)
        self.ax.set_aspect('equal')
        circle1 = plt.Circle((self.x1, self.y1), self.r1, color='red', fill=True, alpha=0.5)
        circle2 = plt.Circle((self.x2, self.y2), self.r2, color='blue', fill=True, alpha=0.5)
        self.ax.add_patch(circle1)
        self.ax.add_patch(circle2)
        self.body1 = circle1
        self.body2 = circle2
        self.trail1, = self.ax.plot([self.x1], [self.y1], 'r-', alpha=0.5)
        self.trail2, = self.ax.plot([self.x2], [self.y2], 'b-', alpha=0.5)

        # Рассчитать и сохранить начальную полную энергию
        self.initial_energy = self.calculate_total_energy()

        self.timer.start(10)

    def update_simulation(self):
        try:
            self.dt = 0.01

            self.x1 += self.vx1 * self.dt
            self.y1 += self.vy1 * self.dt

            self.x2 += self.vx2 * self.dt
            self.y2 += self.vy2 * self.dt

            self.check_wall_collision()
            self.check_body_collision()

            self.body1.center = (self.x1, self.y1)
            self.body2.center = (self.x2, self.y2)

            xdata1, ydata1 = self.trail1.get_data()
            xdata1 = np.append(xdata1, self.x1)
            ydata1 = np.append(ydata1, self.y1)
            self.trail1.set_data(xdata1, ydata1)

            xdata2, ydata2 = self.trail2.get_data()
            xdata2 = np.append(xdata2, self.x2)
            ydata2 = np.append(ydata2, self.y2)
            self.trail2.set_data(xdata2, ydata2)

            self.canvas.draw()

            # Проверка сохранения энергии
            self.check_energy_conservation()

        except RuntimeError:
            self.show_error("Симуляция была прервана.")

    def check_wall_collision(self):
        if self.x1 - self.r1 <= 0:
            self.vx1 = abs(self.vx1)
        elif self.x1 + self.r1 >= self.width:
            self.vx1 = -abs(self.vx1)
        if self.y1 - self.r1 <= 0:
            self.vy1 = abs(self.vy1)
        elif self.y1 + self.r1 >= self.height:
            self.vy1 = -abs(self.vy1)

        if self.x2 - self.r2 <= 0:
            self.vx2 = abs(self.vx2)
        elif self.x2 + self.r2 >= self.width:
            self.vx2 = -abs(self.vx2)
        if self.y2 - self.r2 <= 0:
            self.vy2 = abs(self.vy2)
        elif self.y2 + self.r2 >= self.height:
            self.vy2 = -abs(self.vy2)

    def check_body_collision(self):
        dist = np.sqrt((self.x2 - self.x1) ** 2 + (self.y2 - self.y1) ** 2)
        if dist <= self.r1 + self.r2:
            normal = np.array([self.x2 - self.x1, self.y2 - self.y1]) / dist
            tangent = np.array([-normal[1], normal[0]])

            v1 = np.array([self.vx1, self.vy1])
            v2 = np.array([self.vx2, self.vy2])
            v1n = np.dot(v1, normal)
            v1t = np.dot(v1, tangent)
            v2n = np.dot(v2, normal)
            v2t = np.dot(v2, tangent)

            m1, m2 = self.m1, self.m2
            v1n_new = (v1n * (m1 - m2) + 2 * m2 * v2n) / (m1 + m2)
            v2n_new = (v2n * (m2 - m1) + 2 * m1 * v1n) / (m1 + m2)

            v1n_vec = v1n_new * normal
            v1t_vec = v1t * tangent
            v2n_vec = v2n_new * normal
            v2t_vec = v2t * tangent

            v1_new = v1n_vec + v1t_vec
            v2_new = v2n_vec + v2t_vec
            self.vx1, self.vy1 = v1_new
            self.vx2, self.vy2 = v2_new

            overlap = self.r1 + self.r2 - dist
            correction_vector = normal * (overlap / 2)
            self.x1 -= correction_vector[0]
            self.y1 -= correction_vector[1]
            self.x2 += correction_vector[0]
            self.y2 += correction_vector[1]

    def calculate_total_energy(self):
        # Рассчитать полную кинетическую энергию системы
        kinetic_energy1 = 0.5 * self.m1 * (self.vx1 ** 2 + self.vy1 ** 2)
        kinetic_energy2 = 0.5 * self.m2 * (self.vx2 ** 2 + self.vy2 ** 2)
        total_energy = kinetic_energy1 + kinetic_energy2
        return total_energy

    def check_energy_conservation(self):
        # Проверка сохранения полной кинетической энергии
        current_energy = self.calculate_total_energy()
        energy_loss = abs(current_energy - self.initial_energy) / self.initial_energy

        # Если потеря энергии больше 1%, выводим предупреждение
        if energy_loss > 0.01:
            self.show_error(f'Потеря энергии превышает допустимый порог: {energy_loss:.2%}')

    def closeEvent(self, event):
        self.timer.stop()
        event.accept()


def main():
    app = QApplication(sys.argv)
    window = ElasticCollisionSimulation()
    window.show()
    try:
        sys.exit(app.exec_())
    except SystemExit:
        pass


if __name__ == '__main__':
    main()
