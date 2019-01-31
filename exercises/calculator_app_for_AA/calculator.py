# coding: utf-8

import math
import kivy
kivy.require("1.10.1")
from kivy.app import App
from kivy.uix.boxlayout import BoxLayout
from kivy.properties import BooleanProperty, ListProperty
from kivy.core.text import LabelBase, DEFAULT_FONT
from kivy.resources import resource_add_path
from kivy.core.window import Window
from kivy.utils import get_color_from_hex
from kivy.uix.carousel import Carousel

from kivy.resources import resource_add_path
import sys
if hasattr(sys, "_MEIPASS"):
    resource_add_path(sys._MEIPASS)


class Calculator1(BoxLayout):
    def __init__(self, **kwargs):
        super(Calculator1, self).__init__(**kwargs)

    clear_bool = BooleanProperty(False)

    def print_number(self, number):
        if self.clear_bool:
            self.clear()

        text = "{}{}".format(self.ids["display_input1"].text, number) # ���܂ł̓��͂��ꂽ������Ɠ��͂��ꂽ�l���\�������
        self.ids["display_input1"].text = text

        print("{0} pushed".format(number))

    def print_operator(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        text = "{} {} ".format(self.ids["display_input1"].text, operator)
        self.ids["display_input1"].text = text

        print("{0} pushed".format(operator))

    def change_encode(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        if self.ids["display_input1"].text[:1] == "-":
            text = self.ids["display_input1"].text[1:]
        elif self.ids["display_input1"].text[:2] == " -":
            text = self.ids["display_input1"].text[2:]
        else:
            text = "-" + self.ids["display_input1"].text

        self.ids["display_input1"].text = text
        print("{0} pushed".format(operator))

    def copy_paste(self):
        if self.ids["display_input1"].text == "":
            self.ids["display_input1"].text = CalculatorRoot.copy
        else:
            CalculatorRoot.copy = self.ids["display_input1"].text

        print("C/P pushed")

    def clear(self):
        self.ids["display_input1"].text = ""
        self.clear_bool = False
        print("C pushed")

    def delete(self):
        self.ids["display_input1"].text = self.ids["display_input1"].text[:-1]
        print("DEL pushed")

    def calculate(self):
        try:
            ans=eval(self.ids["display_input1"].text)
            if isinstance(ans, float):
                self.ids["display_input1"].text = "{0:.8f}".format(ans)
            else:
                self.ids["display_input1"].text = str(ans)
            self.clear_bool = True
            print('Cal Completed')
        except:
            print('Input Error')




class Calculator2(BoxLayout):
    def __init__(self, **kwargs):
        super(Calculator2, self).__init__(**kwargs)

    clear_bool = BooleanProperty(False)
    formula = ""
    second_value = 0

    def print_number(self, number):
        if self.clear_bool:
            self.clear()

        self.formula = "{}{}".format(self.formula, number)
        text = "{}{}".format(self.ids["display_input2"].text, number) # ���܂ł̓��͂��ꂽ������Ɠ��͂��ꂽ�l���\�������
        self.ids["display_input2"].text = text
        print("{0} pushed".format(number))

    def print_number_pi(self, number):
        if self.clear_bool:
            self.clear()

        self.formula = "{}{}".format(self.formula, math.pi)
        text = "{}{}".format(self.ids["display_input2"].text, number) # ���܂ł̓��͂��ꂽ������Ɠ��͂��ꂽ�l���\�������
        self.ids["display_input2"].text = text
        print("{0} pushed".format(number))

    def print_operator(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = "{}{}".format(self.formula, operator)
        text = "{} {} ".format(self.ids["display_input2"].text, operator)
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_factorial(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.factorial(eval(self.ids["display_input2"].text)))
        text = "{} {} ".format(self.ids["display_input2"].text, "!")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_root(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.sqrt(eval(self.ids["display_input2"].text)))
        text = "{} {} ".format(operator, self.ids["display_input2"].text)
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_sq(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.pow(eval(self.ids["display_input2"].text), 2))
        text = "{} {} ".format(self.ids["display_input2"].text, operator)
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_pow(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.first_value = eval(self.ids["display_input2"].text)
        print("{0} pushed".format(operator))
        self.all_clear()
        self.second_value = 1

    def print_operator_exp(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.exp(eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_ten_exp(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.pow(10, eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_log_ten(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.log10(eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_ln(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.log(eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_sin(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.sin(eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_cos(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.cos(eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def print_operator_tan(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(math.tan(eval(self.ids["display_input2"].text)))
        text = "{} {} {} {} ".format(operator, "(", self.ids["display_input2"].text, ")")
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def reverse(self):
        if self.clear_bool:
            self.clear_bool = False

        self.formula = str(1/float((self.ids["display_input2"].text)))
        self.ids["display_input2"].text = self.formula
        print("1/x pushed")

    def change_encode(self, operator):
        if self.clear_bool:
            self.clear_bool = False

        if self.ids["display_input2"].text[:1] == "-":
            text = self.ids["display_input2"].text[1:]
        elif self.ids["display_input2"].text[:2] == " -":
            text = self.ids["display_input2"].text[2:]
        else:
            text = "-" + self.ids["display_input2"].text

        self.formula = text
        self.ids["display_input2"].text = text
        print("{0} pushed".format(operator))

    def change_deg(self, operator):
        if self.clear_bool:
            self.clear_bool = False
        deg = eval(self.formula)
        """
        if isinstance(deg, float): #rad
            self.ids["display_input2"].text = math.degrees(deg)
        else: #deg
            self.ids["display_input2"].text = math.radians(deg)
        """
        if deg <10:
            self.formula = str(math.degrees(deg))
            self.ids["display_input2"].text = str(math.degrees(deg))
        else:
            self.formula = str(math.radians(deg))
            self.ids["display_input2"].text = str(math.radians(deg))
        print("{0} pushed".format(operator))

    def copy_paste(self):
        if self.ids["display_input2"].text == "":
            self.ids["display_input2"].text = CalculatorRoot.copy
        else:
            CalculatorRoot.copy = self.ids["display_input2"].text
        print("C/P pushed")

    def clear(self):
        self.ids["display_input2"].text = ""
        self.formula = ""
        self.clear_bool = False
        print("C pushed")

    def delete(self):
        self.ids["display_input2"].text = self.ids["display_input2"].text[:-1]
        self.formula = self.formula[:-1]
        print("DEL pushed")

    def calculate(self):
        if self.second_value == 0:
            try:
                ans=eval(self.formula)
                if isinstance(ans, float):
                    self.ids["display_input2"].text = "{0:.8f}".format(ans)
                else:
                    self.ids["display_input2"].text = str(ans)
                self.clear_bool = True
                print('Cal Completed')
            except:
                print('Input Error')
        else:
            try:
                print("2val ver")
                ans = math.pow(self.first_value, eval(self.formula))
                if isinstance(ans, float):
                    self.ids["display_input2"].text = "{0:.8f}".format(ans)
                else:
                    self.ids["display_input2"].text = str(ans)
                self.clear_bool = True
                self.second_value = 0
                print('Cal Completed')
            except:
                print('Input Error')
            



class Calculator3(BoxLayout):
    def __init__(self, **kwargs):
        super(Calculator3, self).__init__(**kwargs)

    display = "display_input30"
    flag=0

    color_rocket = ListProperty([0,1,0,1])
    color_plane = ListProperty([1,1,1,1])
    color_Vj = ListProperty([1,1,1,1])
    color_Isp = ListProperty([1,1,1,1])
    color_payload_ratio = ListProperty([1,1,1,1])
    color_Cl = ListProperty([1,1,1,1])
    color_Cd = ListProperty([1,1,1,1])
    color_Tn = ListProperty([1,1,1,1])
    color_Rmax = ListProperty([1,1,1,1])

    color30 = ListProperty([1,1,0.7,1])
    color310 = ListProperty([1,1,1,1])
    color320 = ListProperty([1,1,1,1])
    color330 = ListProperty([1,1,1,1])
    color340 = ListProperty([1,1,1,1])
    color350 = ListProperty([1,1,1,1])
    color360 = ListProperty([1,1,1,1])
    color311 = ListProperty([1,1,1,1])
    color321 = ListProperty([1,1,1,1])
    color331 = ListProperty([1,1,1,1])
    color341 = ListProperty([1,1,1,1])
    color351 = ListProperty([1,1,1,1])
    color361 = ListProperty([1,1,1,1])
    color371 = ListProperty([1,1,1,1])
    color381 = ListProperty([1,1,1,1])
    color391 = ListProperty([1,1,1,1])
    color3101 = ListProperty([1,1,1,1])

    def print_number(self, number):
        text = "{}{}".format(self.ids[self.display].text, number) # ���܂ł̓��͂��ꂽ������Ɠ��͂��ꂽ�l���\�������
        self.ids[self.display].text = text
        print("{0} pushed".format(number))

    def copy_paste(self):
        if self.ids[self.display].text == "":
            self.ids[self.display].text = CalculatorRoot.copy
        else:
            CalculatorRoot.copy = self.ids[self.display].text
        print("C/P pushed")

    def all_clear(self):
        if self.flag == 0:
            self.ids["display_input30"].text = ""
            self.ids["display_input310"].text = ""
            self.ids["display_input320"].text = ""
            self.ids["display_input330"].text = ""
            self.ids["display_input340"].text = ""
            self.ids["display_input350"].text = ""
            self.ids["display_input360"].text = ""
        else:
            self.ids["display_input30"].text = ""
            self.ids["display_input311"].text = ""
            self.ids["display_input321"].text = ""
            self.ids["display_input331"].text = ""
            self.ids["display_input341"].text = ""
            self.ids["display_input351"].text = ""
            self.ids["display_input361"].text = ""
            self.ids["display_input371"].text = ""
            self.ids["display_input381"].text = ""
            self.ids["display_input391"].text = ""
            self.ids["display_input3101"].text = ""
        print("AC pushed")

    def clear(self):
        self.ids[self.display].text = ""
        print("C pushed")

    def delete(self):
        self.ids[self.display].text = self.ids[self.display].text[:-1]
        print("DEL pushed")

    def calc_Vj(self):
        try:
            Vj = math.sqrt(2*float(self.ids["display_input350"].text)*3378.59*float(self.ids["display_input330"].text) * (1- math.pow(float(self.ids["display_input340"].text)/float(self.ids["display_input320"].text), 0.365145/1.365145)))
            self.ids["display_input30"].text = "{0:.1f}".format(Vj)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 0
        self.color_Vj = [1,0.7,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
    def calc_Isp(self):
        try:
            Vj = math.sqrt(2*float(self.ids["display_input350"].text)*3378.59*float(self.ids["display_input330"].text) * (1- math.pow(float(self.ids["display_input340"].text)/float(self.ids["display_input320"].text), 0.365145/1.365145)))
            Isp = Vj/9.80665
            self.ids["display_input30"].text = "{0:.1f}".format(Isp)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 0
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,0.7,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
    def calc_payload_ratio(self):
        try:
            Vj = math.sqrt(2*float(self.ids["display_input350"].text)*3378.59*float(self.ids["display_input330"].text) * (1- math.pow(float(self.ids["display_input340"].text)/float(self.ids["display_input320"].text), 0.365145/1.365145)))
            Isp = Vj/9.80665
            rho = (1+float(self.ids["display_input310"].text))/(1/0.071 + float(self.ids["display_input310"].text)/1.14)
            f = 1/(1+30*rho)
            payload_ratio = (math.exp(-float(self.ids["display_input360"].text)/(9.80665*Isp))-f)/(1-f)
            self.ids["display_input30"].text = "{0:.6f}".format(payload_ratio)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 0
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,0.7,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]

    def calc_Cl(self):
        try:
            Cl = 2*float(self.ids["display_input311"].text)*9.80665 /(float(self.ids["display_input341"].text)*math.pow(float(self.ids["display_input321"].text), 2) * float(self.ids["display_input331"].text))
            self.ids["display_input30"].text = "{0:.6f}".format(Cl)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 1
        self.color_Cl = [1,0.7,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
    def calc_Cd(self):
        try:
            Cl = 2*float(self.ids["display_input311"].text)*9.80665 /(float(self.ids["display_input341"].text)*math.pow(float(self.ids["display_input321"].text), 2) * float(self.ids["display_input331"].text))
            Cd = float(self.ids["display_input351"].text) + pow(Cl,2)/(float(self.ids["display_input361"].text)*math.pi*float(self.ids["display_input371"].text))
            self.ids["display_input30"].text = "{0:.6f}".format(Cd)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 1
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,0.7,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
    def calc_Tn(self):
        try:
            Cl = 2*float(self.ids["display_input311"].text)*9.80665 /(float(self.ids["display_input341"].text) * math.pow(float(self.ids["display_input321"].text), 2) * float(self.ids["display_input331"].text))
            Cd = float(self.ids["display_input351"].text) + pow(Cl,2)/(float(self.ids["display_input361"].text)*math.pi*float(self.ids["display_input371"].text))
            Tn = float(self.ids["display_input341"].text) * math.pow(float(self.ids["display_input321"].text), 2) * float(self.ids["display_input331"].text) * Cd/2
            self.ids["display_input30"].text = "{0:.0f}".format(Tn)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 1
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,0.7,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
    def calc_Rmax(self):
        try:
            Rfj = float(self.ids["display_input381"].text)*3600/float(self.ids["display_input391"].text) * math.pow(2/(float(self.ids["display_input341"].text)*float(self.ids["display_input331"].text)), 1/2) * 3/4*math.pow(float(self.ids["display_input361"].text)*math.pi*float(self.ids["display_input371"].text)/(3*math.pow(float(self.ids["display_input351"].text),3)), 1/4)
            Rmax = 2*Rfj*(math.sqrt(float(self.ids["display_input311"].text)*9.80665) - math.sqrt((float(self.ids["display_input311"].text) - float(self.ids["display_input3101"].text))*9.80665))
            self.ids["display_input30"].text = "{0:.0f}".format(Rmax)
        except:
            self.ids["display_input30"].text = "ValueError"
        self.display30()
        self.flag = 1
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,0.7,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]


    def color_reset_rocket(self):
        self.color310 = [1,1,1,1]
        self.color320 = [1,1,1,1]
        self.color330 = [1,1,1,1]
        self.color340 = [1,1,1,1]
        self.color350 = [1,1,1,1]
        self.color360 = [1,1,1,1]
    def color_reset_plane(self):
        self.color311 = [1,1,1,1]
        self.color321 = [1,1,1,1]
        self.color331 = [1,1,1,1]
        self.color341 = [1,1,1,1]
        self.color351 = [1,1,1,1]
        self.color361 = [1,1,1,1]
        self.color371 = [1,1,1,1]
        self.color381 = [1,1,1,1]
        self.color391 = [1,1,1,1]
        self.color3101 = [1,1,1,1]

    def display30(self):
        self.display = "display_input30"
        self.color30 = [1,1,0.7,1]
        self.color_reset_rocket()
        self.color_reset_plane()

    def display310(self):
        self.display = "display_input310"
        self.flag = 0
        self.color_rocket = [0,1,0,1]
        self.color_plane = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color310 = [1,1,0.7,1]
    def display320(self):
        self.display = "display_input320"
        self.flag = 0
        self.color_rocket = [0,1,0,1]
        self.color_plane = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color320 = [1,1,0.7,1]
    def display330(self):
        self.display = "display_input330"
        self.flag = 0
        self.color_rocket = [0,1,0,1]
        self.color_plane = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color330 = [1,1,0.7,1]
    def display340(self):
        self.display = "display_input340"
        self.flag = 0
        self.color_rocket = [0,1,0,1]
        self.color_plane = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color340 = [1,1,0.7,1]
    def display350(self):
        self.display = "display_input350"
        self.flag = 0
        self.color_rocket = [0,1,0,1]
        self.color_plane = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color350 = [1,1,0.7,1]
    def display360(self):
        self.display = "display_input360"
        self.flag = 0
        self.color_rocket = [0,1,0,1]
        self.color_plane = [1,1,1,1]
        self.color_Cl = [1,1,1,1]
        self.color_Cd = [1,1,1,1]
        self.color_Tn = [1,1,1,1]
        self.color_Rmax = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color360 = [1,1,0.7,1]
        
    def display311(self):
        self.display = "display_input311"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color311 = [1,1,0.7,1]
    def display321(self):
        self.display = "display_input321"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color321 = [1,1,0.7,1]
    def display331(self):
        self.display = "display_input331"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color331 = [1,1,0.7,1]
    def display341(self):
        self.display = "display_input341"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color341 = [1,1,0.7,1]
    def display351(self):
        self.display = "display_input351"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color351 = [1,1,0.7,1]
    def display361(self):
        self.display = "display_input361"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color361 = [1,1,0.7,1]
    def display371(self):
        self.display = "display_input371"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color371 = [1,1,0.7,1]
    def display381(self):
        self.display = "display_input381"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color381 = [1,1,0.7,1]
    def display391(self):
        self.display = "display_input391"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color391 = [1,1,0.7,1]
    def display3101(self):
        self.display = "display_input3101"
        self.flag = 1
        self.color_rocket = [1,1,1,1]
        self.color_plane = [0,1,0,1]
        self.color_Vj = [1,1,1,1]
        self.color_Isp = [1,1,1,1]
        self.color_payload_ratio = [1,1,1,1]
        self.color30 = [1,1,1,1]
        self.color_reset_rocket()
        self.color_reset_plane()
        self.color3101 = [1,1,0.7,1]




class Calculator4(BoxLayout):
    def __init__(self, **kwargs):
        super(Calculator4, self).__init__(**kwargs)

    display = "display_input41"
    color41 = ListProperty([1,1,0.7,1])
    color42 = ListProperty([1,1,1,1])
    operator = 0


    def print_number(self, number):
        text = "{}{}".format(self.ids[self.display].text, number) # ���܂ł̓��͂��ꂽ������Ɠ��͂��ꂽ�l���\�������
        self.ids[self.display].text = text
        print("{0} pushed".format(number))
    def copy_paste(self):
        if self.ids[self.display].text == "":
            self.ids[self.display].text = CalculatorRoot.copy
        else:
            CalculatorRoot.copy = self.ids[self.display].text
        print("C/P pushed")
    def all_clear(self):
        self.ids["display_input41"].text = ""
        self.ids["display_input42"].text = ""
        print("AC pushed")
    def clear(self):
        self.ids[self.display].text = ""
        print("C pushed")
    def delete(self):
        self.ids[self.display].text = self.ids[self.display].text[:-1]
        print("DEL pushed")

    def display41(self):
        self.display = "display_input41"
        self.color41 = [1,1,0.7,1]
        self.color42 = [1,1,1,1]
    def display42(self):
        self.display = "display_input42"
        self.color41 = [1,1,1,1]
        self.color42 = [1,1,0.7,1]
    def convert(self):
        if self.display == "display_input41":
            ans = float(self.ids["display_input41"].text) * self.operator
            self.ids["display_input42"].text = "{0:.8f}".format(ans)
            self.display42()
        else:
            ans = float(self.ids["display_input42"].text) / self.operator
            self.ids["display_input41"].text = "{0:.8f}".format(ans)
            self.display41()

    def switch(self):
        a = self.ids["display_input41"].text
        self.ids["display_input41"].text = self.ids["display_input42"].text
        self.ids["display_input42"].text = a
        if self.display == "display_input41":
            self.display42()
        else:
            self.display41()

    def camma_up41(self):
        self.ids["display_input41"].text = str(1000*(float(self.ids["display_input41"].text)))
    def camma_down41(self):
        self.ids["display_input41"].text = str(0.001*(float(self.ids["display_input41"].text)))
    def camma_up42(self):
        self.ids["display_input42"].text = str(1000*(float(self.ids["display_input42"].text)))
    def camma_down42(self):
        self.ids["display_input42"].text = str(0.001*(float(self.ids["display_input42"].text)))

    def ft_m(self):
        self.operator = 0.3048
        self.ids["unit1"].text = "[ft]"
        self.ids["unit2"].text = "[m]"
    def nm_km(self):
        self.operator = 1.852
        self.ids["unit1"].text = "[nm]"
        self.ids["unit2"].text = "[km]"
    def kt_kmh(self):
        self.operator = 1.852
        self.ids["unit1"].text = "[kt]"
        self.ids["unit2"].text = "[km/h]"
    def ms_kmh(self):
        self.operator = 3.6
        self.ids["unit1"].text = "[m/s]"
        self.ids["unit2"].text = "[km/h]"
    def gallon_L(self):
        self.operator = 3.7854
        self.ids["unit1"].text = "[gal]"
        self.ids["unit2"].text = "[L]"
    def lb_kg(self):
        self.operator = 0.454
        self.ids["unit1"].text = "[lb]"
        self.ids["unit2"].text = "[kg]"
    def lb_slug(self):
        self.operator = 0.03108486
        self.ids["unit1"].text = "[lb]"
        self.ids["unit2"].text = "[slug]"
    def atm_Pa(self):
        self.operator = 101300
        self.ids["unit1"].text = "[atm]"
        self.ids["unit2"].text = "[Pa]"
    def psi_Pa(self):
        self.operator = 6894
        self.ids["unit1"].text = "[psi]"
        self.ids["unit2"].text = "[Pa]"
    def HP_kW(self):
        self.operator = 0.7457
        self.ids["unit1"].text = "[HP]"
        self.ids["unit2"].text = "[kW]"
    def PS_kW(self):
        self.operator = 0.7355
        self.ids["unit1"].text = "[PS]"
        self.ids["unit2"].text = "[kW]"






class CalculatorRoot(BoxLayout):
    color1 = ListProperty([0,1,0,1])
    color2 = ListProperty([1,1,1,1])
    color3 = ListProperty([1,1,1,1])
    color4 = ListProperty([1,1,1,1])

    copy = ""

    def __init__(self, **kwargs):
        super(CalculatorRoot, self).__init__(**kwargs)

    def change_color1(self):
        self.color1 = [0,1,0,1]
        self.color2 = [1,1,1,1]
        self.color3 = [1,1,1,1]
        self.color4 = [1,1,1,1]
    def change_color2(self):
        self.color1 = [1,1,1,1]
        self.color2 = [0,1,0,1]
        self.color3 = [1,1,1,1]
        self.color4 = [1,1,1,1]
    def change_color3(self):
        self.color1 = [1,1,1,1]
        self.color2 = [1,1,1,1]
        self.color3 = [0,1,0,1]
        self.color4 = [1,1,1,1]
    def change_color4(self):
        self.color1 = [1,1,1,1]
        self.color2 = [1,1,1,1]
        self.color3 = [1,1,1,1]
        self.color4 = [0,1,0,1]

    def change_calc1(self):   
        self.carousel.load_slide(self.calculator1)
        self.change_color1()
        self.current_mode = 1
    def change_calc2(self):   
        self.carousel.load_slide(self.calculator2)
        self.change_color2()
        self.current_mode = 2
    def change_calc3(self):   
        self.carousel.load_slide(self.calculator3)
        self.change_color3()
        self.current_mode = 3
    def change_calc4(self):   
        self.carousel.load_slide(self.calculator4)
        self.change_color4()
        self.current_mode = 4

    def on_index(self, instance, value): #�X���C�v�Ń��[�h�ς�����Ƃ��ɐF�\�����t���čs������
        cur_mode = value+1
        print("current mode is", cur_mode)
        if cur_mode == 1:
            self.change_color1()
        elif cur_mode == 2:
            self.change_color2()
        elif cur_mode == 3:
            self.change_color3()
        else:
            self.change_color4()



class CalculatorApp(App):
    def __init__(self, **kwargs):
        super(CalculatorApp, self).__init__(**kwargs)
        self.title = "Calculator"
    pass

if __name__ == "__main__":
    CalculatorApp().run()