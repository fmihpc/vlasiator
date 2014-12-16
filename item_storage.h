/*
 * This file is part of Vlasiator.
 * Copyright 2013, 2014 Finnish Meteorological Institute
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
