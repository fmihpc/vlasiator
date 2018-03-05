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

#ifndef OBJECT_FACTORY_H
#define OBJECT_FACTORY_H

#include <iostream>
#include <map>

#include "definitions.h"

/** A generic storage class for storing items, such as variable values or 
 * function pointers. The storage class can store any type of item that does 
 * not require a constructor call, i.e., it does not create a _new_ copy of 
 * an item when one is requested.*/
template<typename PRODUCT>
class ObjectFactory {
 public:
   
   PRODUCT* create(const std::string& name) const;
   bool add(const std::string& name,PRODUCT* (*maker)());
   size_t size() const;

 private:

   // Here the mysterious "PRODUCT* (*)()" is just a function pointer.
   // If it had a name 'maker', it could be expanded as
   // PRODUCT* (*maker)()
   // In other words, it is a pointer to a function that takes no arguments, 
   // and that returns a pointer to PRODUCT.
   std::map<std::string,PRODUCT* (*)() > manufacturers; /**< Container for all stored maker functions.*/
};

/** Get an item (product) with the given name.
 * @param name The name of the item.
 * @param func The requested item is stored here, if it was found.
 * @return If true, variable func contains a valid item.*/
template<typename PRODUCT> inline
PRODUCT* ObjectFactory<PRODUCT>::create(const std::string& name) const {
   typename std::map<std::string,PRODUCT* (*)()>::const_iterator it = manufacturers.find(name);
   if (it == manufacturers.end()) return NULL;
   return (*it->second)();
}

/** Register a maker to the factory. This function will fail 
 * to succeed if the factory already contains a maker with the given name.
 * @param name A unique name for the maker.
 * @param func The maker function.
 * @return If true, the maker was added to the factory.*/
template<typename PRODUCT> inline
bool ObjectFactory<PRODUCT>::add(const std::string& name,PRODUCT* (*maker)()) {
   // The insert returns a pair<iterator,bool>, where the boolean value is 'true' 
   // if the maker function was inserted to the map. Skip the pair creation 
   // and just return the boolean.
   return manufacturers.insert(make_pair(name,maker)).second;
}

template<typename PRODUCT> inline
size_t ObjectFactory<PRODUCT>::size() const {
   return this->manufacturers.size();
}

#endif
