# Implementations of geomagnetic storm list algorithms.
# Copyright (C) 2025 John Coxon (work@johncoxon.co.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from .list import List
from pandas import DataFrame

class WG19(List):
    """
    Geomagnetic storms identified using the algorithm described in Walach and Grocott (2019).

    ----------
    category_bins : np.ndarray[int]
    end_time : datetime.datetime
        The end of the window in which storms were identified.
    quiet_level : int
    start_time : datetime.datetime
        The start of the window in which storms were identified.
    list : pandas.DataFrame
        A DataFrame describing each of the identified storms.
    """

    def __init__(self, time, symh, quiet_level=-15, category_bins=np.array([-80, -150, -300])):
        """
        Parameters
        ----------
        time : np.ndarray[datetime.datetime]
            Timestamps of the OMNI data being used.
        symh : np.ndarray[int]
            SYM-H levels in nT from OMNI.
        quiet_level : int, optional, default -15
            The quiet level to use, defaults to -15 nT as in the original paper.
        category_bins : np.ndarray[int], optional, default [-50, -150, -300]
            The thresholds to use for the three storm categories, defaults to the thresholds used in the original paper.
        """
        self.category_bins = category_bins
        self.quiet_level = quiet_level

        self._crossings = self._find_quiet_level_crossings(symh)
        self._storm_categories = np.zeros_like(time, dtype=int)

        super().__init__(time, self.find_storms(time, symh))

    def _find_quiet_level_crossings(self, symh):
        """
        Find crossings in the SYM-H index by subtracting the quiet level from the data and finding zero crossings.

        All points equal to the quiet level are identified as crossings, and all points BEFORE and AFTER a crossing
        through the quiet level are identified as crossings.

        For example, if quiet level is set to -15 nT and the data goes from -17 to -13 nT, the timestamps associated
        with both the -17 nT value and the -13 nT value will be identified as crossing points.

        Parameters
        ----------
        symh : np.ndarray[int]

        Returns
        -------
        crossings : np.ndarray[int]
            An array of integers describing the points at which SYM-H goes through the quiet level.
        """
        through_zero = np.where(symh == self.quiet_level)[0]

        # Find turning points which don't go through the quiet level, so e.g. from 1 to -1 without passing through 0.
        # This line identifies the points immediately before the zero crossing.
        modified_symh = symh - self.quiet_level
        crossings_left = np.where(np.diff(np.sign(modified_symh)))[0]

        # Add 1 to the above array so that we flag the timestamps BOTH before AND after the crossing as turning points.
        # (There will be an error of <1 minute no matter which way we resolve the ambiguity and this way is easier.)
        crossings_right = crossings_left + 1

        crossings = np.sort(np.concatenate((through_zero, crossings_left, crossings_right)))

        return crossings

    def find_storms(self, time, symh):
        """
        Find all the storms in the given timeframe.

        Parameters
        ----------
        time : np.ndarray[datetime.datetime]
        symh : np.ndarray[int]

        Returns
        -------
        storms : np.ndarray[Storm]
            The list of storms in the given timeframe.
        """
        storms = []

        while True:
            storms.append(self.find_storm(time, symh))
            if np.nanmin(symh[self._storm_categories == 0]) >= np.max(self.category_bins):
                break

        storms = DataFrame(
            {"initial_phase": np.array([storm["initial_start"] for storm in storms]),
             "main_phase": np.array([storm["main_start"] for storm in storms]),
             "recovery_phase": np.array([storm["recovery_start"] for storm in storms]),
             "end": np.array([storm["end"] for storm in storms]),
             "minimum_symh": np.array([storm["minimum_symh"] for storm in storms]).astype(int),
             "category": np.array([storm["category"] for storm in storms]).astype(int)}
        ).sort_values("recovery_phase").reset_index()
        storms.pop("index")

        return storms

    def find_storm(self, time, symh):
        """
        Find a storm that has not yet been identified.

        Parameters
        ----------
        time : np.ndarray[datetime.datetime]
        symh : np.ndarray[int]

        Returns
        -------
        storm : Storm
            The next identified storm.
        """
        minimum_symh = np.nanmin(symh[self._storm_categories == 0])

        recovery_start = np.where((symh == minimum_symh) & (self._storm_categories == 0))[0][0].astype(int)
        storm_category = np.digitize(minimum_symh, self.category_bins)

        end_turning_point = np.searchsorted(self._crossings, recovery_start)
        main_turning_point = end_turning_point - 1
        main_start = self._crossings[main_turning_point]

        storm_end = self._crossings[end_turning_point]
        if storm_end == recovery_start:
            storm_end = self._crossings[end_turning_point + 1]
        if storm_end >= len(time):
            storm_end = len(time) - 1

        n_minutes = 18 * 60
        initial_peak = (np.nanargmax(symh[main_start - n_minutes:main_start])
                        + main_start - n_minutes)

        initial_turning_point = np.searchsorted(self._crossings, initial_peak) - 1
        if initial_turning_point < 0:
            initial_start = 0
        else:
            initial_start = self._crossings[initial_turning_point]

        self._storm_categories[initial_start:storm_end + 1] = storm_category

        storm = {"initial_start": time[initial_start],
                 "main_start": time[main_start],
                 "recovery_start": time[recovery_start],
                 "end": time[storm_end],
                 "minimum_symh": minimum_symh,
                 "category": storm_category}

        return storm

    def _string(self):
        if self.start_time.year == self.end_time.year:
            date_string = f"{self.start_time:%Y}"
        else:
            date_string = f"{self.start_time:%Y} to {self.end_time:%Y}"

        return (
            f"Geomagnetic storms for {date_string} from Walach & Grocott (2019).\n"
            f"Quiet level {self.quiet_level} nT.\nCategory thresholds "
            f"[{self.category_bins[0]} nT, {self.category_bins[1]} nT, {self.category_bins[2]} nT].")

    def __str__(self):
        return self._string()

    def __repr__(self):
        return self._string()
