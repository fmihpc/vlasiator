/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifndef ITEM_STORAGE_H
#define ITEM_STORAGE_H

#include <map>

#include "definitions.h"

/** A generic storage class for storing items, such as variable values or 
 * function pointers. The storage class can store any type of item that does 
 * not require a constructor call, i.e., it does not create a _new_ copy of 
 * an item when one is requested.*/
template<typename ITEM>
class ItemStorage {
 public:
   bool get(const std::string& name,ITEM& item) const;
   bool store(const std::string& name,ITEM item);
   
 private:
   std::map<std::string,ITEM> warehouse; /**< Container for all stored items.*/
};

/** Get an item (product) with the given name.
 * @param name The name of the item.
 * @param func The requested item is stored here, if it was found.
 * @return If true, variable func contains a valid item.*/
template<typename ITEM> inline
bool ItemStorage<ITEM>::get(const std::string& name,ITEM& item) const {
   typename std::map<std::string,ITEM>::const_iterator it = warehouse.find(name);
   if (it == warehouse.end()) return false;
   item = it->second;
   return true;
}

/** Register a item (product) to the factory. This function will fail 
 * to succeed if the factory already contains an item (product) with the same name.
 * @param name A unique name for the item (product).
 * @param func Registered item (product).
 * @return If true, the item was added to the factory.*/
template<typename ITEM> inline
bool ItemStorage<ITEM>::store(const std::string& name,ITEM item) {
   // The insert returns a pair<iterator,bool>, where the boolean value is 'true' 
   // if the function was inserted to the map. Skip the pair creation and just return the boolean.
   return warehouse.insert(make_pair(name,item)).second;
}

#endif
