from urllib.request import urlopen
import pandas as pd
import numpy as np


class getXsections:
    """getXsections

    Downloads cross section information from https://clas.sinp.msu.ru/cgi-bin/jlab/db.cgi and parses into dataframes

    Only tested for getting cross sections currently. Not fulling working API for the website.
    """

    def __init__(self, url=None):
        self._url = url
        self._text = self.pullFromUrl()
        # self.datasets = parseText(self._text)

    def url(self, url=None):
        self._url = url

    def pullFromUrl(self):
        if self._url is not None and self._url != "":
            html = urlopen(self._url).read().decode("utf-8")
        else:
            html = None
        return html

    def parseDataSets(self):
        start_lines = []
        lines = self._text.split("\n")
        for n, line in enumerate(lines):
            if "Measurement" in line:
                start_lines.append(n)
        start_lines.append(len(lines))

        self._datasets = []
        for s, start in enumerate(start_lines[:-1]):
            end = start_lines[s + 1] - 1
            self._datasets.append(dataset(lines[start:end]))

    @property
    def text(self):
        return self._text

    @property
    def datasets(self):
        return self._datasets


class dataset:
    """
    Class with 
    """

    def __init__(self, text):
        super().__init__()
        # Full text of the dataset
        self._text = text
        # Dict of info from the headers
        self._experiment = {
            "Measurement": None,
            "Title": None,
            "Spokespersons": None,
            "Year": None,
            "Q^2 range": None,
            "Q2_min": None,
            "Q2_max": None,
            "W range": None,
            "W_min": None,
            "W_max": None,
            "E_beam": None,
            "coloumnNamesWithUnits": None,
            "coloumnNames": None,
        }

        self.dataFrame = None
        self.parseText()

    def parseText(self):
        # Parse through header to fill experiment info
        for key, value in self._experiment.items():
            # Fill dict based on header info
            for i, line in enumerate(self._text[0:20]):
                if "Data" in line:
                    data_line = i
                if key in line:
                    key_location = int(line.find(key) + len(key) + 1)
                    self._experiment[key] = line[key_location:].strip("\r")

        # Get numberical values from w range by finding dash and getting np.float before and after
        if self._experiment["W range"] is not None:
            minus = self._experiment["W range"].find("-")
            self._experiment["W_min"] = np.float(self._experiment["W range"][:minus])
            self._experiment["W_max"] = np.float(
                self._experiment["W range"][minus + 1 :]
            )

        # Get numberical values from w range by finding dash and getting np.float before and after
        if self._experiment["W range"] is not None:
            minus = self._experiment["Q^2 range"].find("-")
            self._experiment["Q2_min"] = np.float(self._experiment["Q^2 range"][:minus])
            self._experiment["Q2_max"] = np.float(
                self._experiment["Q^2 range"][minus + 1 :]
            )

        # Cleanup names and units rows and split to make arrays
        names = self._text[data_line].strip("\r").strip("\n").split("\t")
        units = self._text[data_line + 1].strip("\r").strip("\n").split("\t")
        coloumnNamesWithUnits = []
        coloumnNames = []
        for n, u in zip(names, units):
            coloumnNamesWithUnits.append(f"{n} [{u}]")
            coloumnNames.append(f"{n}")
        self._experiment["coloumnNamesWithUnits"] = coloumnNamesWithUnits
        self._experiment["coloumnNames"] = coloumnNames
        all_data = []
        for data in self._text[data_line + 2 :]:
            split_data = data.strip("\r").split("\t")
            if "" in split_data:
                replace = split_data.index("")
                split_data[replace] = "nan"
            all_data.append(np.array(split_data).astype(np.float))

        self.dataFrame = pd.DataFrame(all_data, columns=self.experiment["coloumnNames"])

    @property
    def experiment(self):
        return self._experiment
