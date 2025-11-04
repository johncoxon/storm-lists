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

class P13(List):
    """
    Geomagnetic storms identified using the algorithm described in Partamies et al. (2013).

    Variables
    ---------
    list : pandas.DataFrame
        A DataFrame describing each of the identified storms.
    start_time : datetime.datetime
        The start of the window in which storms were identified.
    end_time : datetime.datetime
        The end of the window in which storms were identified.
    quiet_level : int
        The quiet level used to generate the storm list.
    category_bins : np.ndarray[int]
        The thresholds used for the three storm categories.
    """

    def __init__(self, time, dst, start_level=0, end_level=-15, minimum_dst=-15, rate_threshold=-2):
        """
        Parameters
        ----------
        time : np.ndarray[datetime.datetime]
            Timestamps of the OMNI data being used.
        dst : np.ndarray[int]
            Dst levels in nT.
        self.start_level : int, optional, default 0
            The start level to use, defaults to 0 nT as in the original paper.
        self.end_level : int, optional, default -15
            The end level to use, defaults to -15 nT as in the original paper.
        minimum_dst : int, optional, default -15
            The minimum level to use, defaults to -15 nT as in the original paper.
        rate_threshold : int, optional, default 2
            The threshold value of dDst/dt in nT per hour, defaults to -2 nT per hour as in the original paper.
        """
        self.start_level = start_level
        self.end_level = end_level
        self.minimum_dst_threshold = minimum_dst
        self.rate_threshold = rate_threshold

        self._rate_of_change = np.concat((np.array([0]), np.diff(dst))).astype(int)

        self._start_crossings = self.find_below_to_above(dst, self.start_level)
        self._end_crossings = self.find_below_to_above(dst, self.end_level)
        self._rate_crossings = self.find_above_to_below(self._rate_of_change, self.rate_threshold)
        self._storm_idents = np.zeros_like(time, dtype=int)

        super().__init__(time, self.find_storms(time, dst))


    @staticmethod
    def find_below_to_above(variable, level):
        """Find times at which the variable goes from under the threshold value to over the threshold value."""
        modified_variable = variable - level

        # Identify the times where the variable is equal to the level, and then only take the ones where it
        # exceeds the level in the subsequent bin.
        crossings_exact = []
        through_zero = np.where(modified_variable == 0)[0]
        for i in through_zero:
            if modified_variable[i + 1] > 0:
                crossings_exact.append(i)

        crossings_exact = np.array(crossings_exact)

        # Find turning points which don't go through the threshold, so e.g. from -1 to 1 without passing through 0.
        # This identifies the points immediately after the crossing.
        crossings_right = np.where(np.diff(np.sign(modified_variable)) == 2)[0] + 1

        crossings = np.sort(np.concatenate([crossings_exact, crossings_right]).astype(int))

        return crossings

    @staticmethod
    def find_above_to_below(variable, level):
        """Find times at which the variable goes from above the threshold value to under the threshold value."""
        modified_variable = variable - level

        # Identify the times where the variable is equal to the level, and then only take the ones where it
        # exceeds the level in the subsequent bin.
        crossings_exact = []
        through_zero = np.where(modified_variable == 0)[0]
        for i in through_zero:
            if modified_variable[i + 1] < 0:
                crossings_exact.append(i)

        crossings_exact = np.array(crossings_exact)

        # Find turning points which don't go through the threshold, so e.g. from -1 to 1 without passing through 0.
        # This identifies the points immediately before the crossing.
        crossings_right = np.where(np.diff(np.sign(modified_variable)) == -2)[0] + 1

        crossings = np.sort(np.concatenate([crossings_exact, crossings_right]).astype(int))

        return crossings

    def find_storms(self, input_time, input_dst):
        """
        Find all the storms in the given input_timeframe.

        Parameters
        ----------
        input_time : np.ndarray[dateinput_time.dateinput_time]
        input_dst : np.ndarray[int]

        Returns
        -------
        storms : pandas.DataFrame
            The list of storms in the given input_timeframe.
        """
        storms = []

        while True:
            next_storm = self.find_storm(input_time, input_dst)
            if next_storm:
                storms.append(next_storm)
            if np.nanmin(input_dst[self._storm_idents == 0]) >= self.minimum_dst_threshold:
                break

        storms = DataFrame(
            {"initial_phase": np.array([storm["initial_start"] for storm in storms]),
             "main_phase": np.array([storm["main_start"] for storm in storms]),
             "recovery_phase": np.array([storm["recovery_start"] for storm in storms]),
             "end": np.array([storm["end"] for storm in storms]),
             "minimum_dst": np.array([storm["minimum_dst"] for storm in storms]).astype(int),
             "storm_type": np.array([storm["storm_type"] for storm in storms]).astype(int)}
        ).sort_values("recovery_phase").reset_index()
        storms.pop("index")

        return storms

    def find_storm(self, input_time, input_dst):
        """
        Find a storm in the data using the Partamies et al. (2013) algorithm.

        Parameters
        ----------
        input_time : np.ndarray[datetime.datetime]
        input_dst : np.ndarray[int]

        Returns
        -------
        storms : dict
            The next identified storm.
        """
        minimum_input_dst = np.nanmin(input_dst[self._storm_idents == 0])
        recovery_start = np.where((input_dst == minimum_input_dst) & (self._storm_idents == 0))[0][0].astype(int)

        # If the rate criterion is violated, flag this minimum as a -1 and bail.
        if self._rate_of_change[recovery_start] >= self.rate_threshold:
            self._storm_idents[recovery_start] = -1
            return None

        rate_crossing = np.searchsorted(self._rate_crossings, recovery_start)

        # If the data ends before a storm is found, it will result in IndexError, so keep an eye out for that.
        if rate_crossing == len(self._rate_crossings):
            self._storm_idents[recovery_start] = -1
            return None

        # This triggers if the hour of minimum Dst is the only one with a rate of change less than 2.
        if self._rate_crossings[rate_crossing] == recovery_start:
            main_start = recovery_start

        # It's possible for the code to find a value of -2, if that happens add one to the main start.
        elif self._rate_of_change[self._rate_crossings[rate_crossing - 1]] == -2:
            main_start = self._rate_crossings[rate_crossing - 1] + 1

        # If that doesn't trigger, everything is kosher.
        else:
            main_start = self._rate_crossings[rate_crossing - 1]

        start_crossing = np.searchsorted(self._start_crossings, main_start) - 1
        initial_start = self._start_crossings[start_crossing]

        # If the initial looping back round, set it to zero.
        if initial_start > main_start:
            initial_start = 0

        # Check whether the start crossing is already identified as a storm. If no, assign as normal.
        if self._storm_idents[initial_start] == 0:
            storm_ident = 1

        # If yes, then do something a little different.
        else:
            not_yet_storm = np.where(self._storm_idents == 0)[0]
            new_initial_start = not_yet_storm[np.searchsorted(not_yet_storm, initial_start)]
            initial_start = new_initial_start
            storm_ident = 2

        # TODO: there is probably an error waiting to happen here where the storm ends after the data ends
        end_crossing = np.searchsorted(self._end_crossings, recovery_start)
        storm_end = self._end_crossings[end_crossing] - 1

        # Do some checks that everything is okay.
        if (self._rate_of_change[main_start:recovery_start] >= self.rate_threshold).any():
            raise ValueError("Entire main phase should have rate of change below the threshold.")
        if self._rate_of_change[main_start - 1] < self.rate_threshold:
            raise ValueError("Rate of change before main phase should be above the threshold.")
        if (input_dst[initial_start] < self.start_level) and (storm_ident == 1) and (initial_start != 0):
            raise ValueError("Start of the storm is below the start threshold.")
        if input_dst[initial_start - 1] > self.start_level:
            if input_dst[initial_start] != 0:
                raise ValueError("Start threshold crossed before the start of the storm.")
        if input_dst[storm_end] > self.end_level:
            raise ValueError(f"End of the storm is above the end threshold.")
        if input_dst[storm_end + 1] < self.end_level:
            raise ValueError(f"End threshold crossed after the end of the storm.")

        self._storm_idents[initial_start:storm_end + 1] = storm_ident

        storm = {"initial_start": input_time[initial_start],
                 "main_start": input_time[main_start],
                 "recovery_start": input_time[recovery_start],
                 "end": input_time[storm_end],
                 "minimum_dst": minimum_input_dst,
                 "storm_type": storm_ident}

        return storm

    def _string(self):
        if self.start_time.year == self.end_time.year:
            date_string = f"{self.start_time:%Y}"
        else:
            date_string = f"{self.start_time:%Y} to {self.end_time:%Y}"

        return (
            f"Geomagnetic storms for {date_string} from Partamies et al. (2013).\n"
            f"Start level {self.start_level} nT.\n"
            f"End level {self.end_level} nT.\n"
            f"Minimum Dst {self.minimum_dst_threshold} nT.\n"
            f"dDst/dt {self.rate_threshold} nT/hour.\n")

    def __str__(self):
        return self._string()

    def __repr__(self):
        return self._string()
