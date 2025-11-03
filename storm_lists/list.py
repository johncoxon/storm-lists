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

class List(object):
    """
    Base class for geomagnetic storm identifications.

    Variables
    ----------
    list : pandas.DataFrame
        A DataFrame describing each of the identified storms.
    start_time : datetime.datetime
        The start of the window in which storms were identified.
    end_time : datetime.datetime
        The end of the window in which storms were identified.
    """
    def __init__(self, time, storm_list):
        """
        Parameters
        ----------
        time : np.ndarray[datetime.datetime]
            Timestamps of the data being used.
        """
        self.list = storm_list
        self.start_time = time[0]
        self.end_time = time[-1]

    def _string(self):
        if self.start_time.year == self.end_time.year:
            date_string = f"{self.start_time:%Y}"
        else:
            date_string = f"{self.start_time:%Y} to {self.end_time:%Y}"

        return f"Geomagnetic storms for {date_string}."

    def __str__(self):
        return self._string()

    def __repr__(self):
        return self._string()
