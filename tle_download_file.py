#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
tle_download_file.py: Script para obter os arquivos .tle a partir do Spacetrack.
"""

import argparse
from datetime import datetime, timedelta
import os
import sys
import time
import configparser
import getpass
from spacetrack import SpaceTrackClient
import spacetrack.operators as op
from alive_progress import alive_bar, config_handler

__author__ = "Filipe Rodrigues"
__version__ = "1.1"
__maintainer__ = "Filipe Rodrigues"
__email__ = "frodriguesfajr@tgmail.com"
__status__ = "Production"


class TleDownload:
    """
    Classe para gerenciar o download de arquivos TLE a partir do Spacetrack.
    """

    DEFAULT_START_DATE = datetime(1957, 10, 4)
    DEFAULT_END_DATE = datetime(2040, 1, 1)
    DEFAULT_NORAD_ID = "56215"
    DEFAULT_DAYS = 2

    def __init__(self):
        self.norad = None
        self.identity = None
        self.password = None
        self.date_start = None
        self.date_end = None
        self.days_number = None
        self.path_tle = None

    @staticmethod
    def parse_args(args=None):
        """
        Parseia os argumentos da linha de comando.
        """
        parser = argparse.ArgumentParser(description="Parses command.")
        parser.add_argument('-credentials_path', help='Caminho para o arquivo de credenciais.')
        parser.add_argument('-identity', help='Seu username do Spacetrack.')
        parser.add_argument('-password', help='Sua senha do Spacetrack.')
        parser.add_argument('-date_start', help='Data inicial no formato DD/MM/AAAA.')
        parser.add_argument('-date_end', help='Data final no formato DD/MM/AAAA.')
        parser.add_argument('-days_number', type=int, help='Número de dias a partir da data atual.')
        parser.add_argument('-norad_id', help='ID NORAD do objeto.')
        return parser.parse_args(args)

    def get_inputs(self, args=None):
        """
        Processa os argumentos e configura as variáveis da classe.
        """
        options = self.parse_args(args)

        # Credenciais
        if options.credentials_path and os.path.isfile(options.credentials_path):
            config = configparser.ConfigParser()
            config.read(options.credentials_path)
            self.identity = config.get("configuration", "username")
            self.password = config.get("configuration", "password")
        else:
            self.identity = options.identity or getpass.getpass('Informe Spacetrack username: ')
            self.password = options.password or getpass.getpass('Informe Spacetrack password: ')

        # Datas
        self.date_start = self.parse_date(options.date_start, self.DEFAULT_START_DATE)
        self.date_end = self.parse_date(options.date_end, self.DEFAULT_END_DATE)

        # Outros parâmetros
        self.days_number = options.days_number or self.DEFAULT_DAYS
        self.norad = options.norad_id or self.DEFAULT_NORAD_ID

        # Caminho para salvar os arquivos
        self.path_tle = os.path.join(os.getcwd(), 'tle_files', self.norad)

    @staticmethod
    def parse_date(date_str, default_date):
        """
        Converte uma “string” de data em objeto datetime.
        """
        if not date_str:
            return default_date
        try:
            return datetime.strptime(date_str, '%d/%m/%Y')
        except ValueError:
            print(f"Formato de data inválido: {date_str}. Usando valor padrão.")
            return default_date

    def search_tle_files(self):
        """
        Realiza o download dos arquivos TLE.
        """
        # Ajusta datas para intervalos padrão, se necessário
        if self.date_end == self.DEFAULT_END_DATE:
            self.date_end = datetime.now()
        if self.date_start == self.DEFAULT_START_DATE:
            self.date_start = datetime.now() - timedelta(days=self.days_number)

        # Cria o diretório, se necessário
        os.makedirs(self.path_tle, exist_ok=True)

        # Conexão com o Spacetrack
        st_client = SpaceTrackClient(self.identity, self.password)
        date_range = op.inclusive_range(self.date_start, self.date_end)
        tle_lines = st_client.tle(
            iter_lines=True,
            norad_cat_id=[int(self.norad)],
            epoch=date_range,
            orderby='epoch',
            format='tle'
        )

        # Processa os arquivos TLE
        tle_data = list(tle_lines)
        for i in range(0, len(tle_data), 2):
            tle_epoch = self.extract_epoch(tle_data[i])
            filename = f"{self.norad}_{tle_epoch}.tle"
            with open(os.path.join(self.path_tle, filename), 'w') as tle_file:
                tle_file.write(tle_data[i] + '\n' + tle_data[i + 1])

    @staticmethod
    def extract_epoch(tle_line):
        """
        Extrai a data de epoch de uma linha TLE.
        """
        year = int(tle_line[18:20])
        year += 2000 if year < 57 else 1900
        day_of_year = float(tle_line[20:32])
        epoch_date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
        return epoch_date.strftime("%Y%m%d_%H%M%S")


# Configurações globais
config_handler.set_global(length=30, spinner='dots_waves')

if __name__ == "__main__":
    tle_downloader = TleDownload()
    with alive_bar(2, bar='smooth') as bar:
        tle_downloader.get_inputs()
        bar()
        tle_downloader.search_tle_files()
        bar()
