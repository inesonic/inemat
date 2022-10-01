##-*-makefile-*-########################################################################################################
# Copyright 2016 - 2022 Inesonic, LLC
# 
# This file is licensed under two licenses.
#
# Inesonic Commercial License, Version 1:
#   All rights reserved.  Inesonic, LLC retains all rights to this software, including the right to relicense the
#   software in source or binary formats under different terms.  Unauthorized use under the terms of this license is
#   strictly prohibited.
#
# GNU Public License, Version 2:
#   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
#   License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
#   version.
#   
#   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
#   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
#   details.
#   
#   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free
#   Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
########################################################################################################################

########################################################################################################################
# Basic build characteristics
#

TEMPLATE = lib

CONFIG -= qt
CONFIG += c++14

win32 {
    # Defined static to make Windows happy.  No binary so it shouldn't matter.  Defaults work for other platforms.
    CONFIG += static
}

unix {
    QMAKE_CFLAGS_DEBUG -= -O
    QMAKE_CFLAGS_DEBUG -= -O1
    QMAKE_CFLAGS_DEBUG -= -O2

    QMAKE_CFLAGS_DEBUG *= -O3

    QMAKE_CXXFLAGS_DEBUG -= -O
    QMAKE_CXXFLAGS_DEBUG -= -O1
    QMAKE_CXXFLAGS_DEBUG -= -O2

    QMAKE_CXXFLAGS_DEBUG *= -O3

    QMAKE_CFLAGS_RELEASE -= -O
    QMAKE_CFLAGS_RELEASE -= -O1
    QMAKE_CFLAGS_RELEASE -= -O2

    QMAKE_CFLAGS_RELEASE *= -O3

    QMAKE_CXXFLAGS_RELEASE -= -O
    QMAKE_CXXFLAGS_RELEASE -= -O1
    QMAKE_CXXFLAGS_RELEASE -= -O2

    QMAKE_CXXFLAGS_RELEASE *= -O3
}

CONFIG(debug, debug|release) {
    DEFINES += DEBUG_BUILD
} else {
    DEFINES += RELEASE_BUILD
}

########################################################################################################################
# Public includes
#

INCLUDEPATH += include
API_HEADERS = include/mat_common.h \
              include/mat_api.h \

########################################################################################################################
# Setup headers and installation
#

HEADERS = $$API_HEADERS

########################################################################################################################
# Locate build intermediate and output products
#

TARGET = inemat

CONFIG(debug, debug|release) {
    unix:DESTDIR = build/debug
    win32:DESTDIR = build/Debug
} else {
    unix:DESTDIR = build/release
    win32:DESTDIR = build/Release
}

OBJECTS_DIR = $${DESTDIR}/objects
MOC_DIR = $${DESTDIR}/moc
RCC_DIR = $${DESTDIR}/rcc
UI_DIR = $${DESTDIR}/ui
