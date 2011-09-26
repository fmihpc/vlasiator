/*
This file is part of Vlasiator.

Copyright 2010, 2011 Finnish Meteorological Institute

Vlasiator is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3
as published by the Free Software Foundation.

Vlasiator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MUXML_H
#define MUXML_H

#include <ostream>
#include <map>
#include <list>
#include <utility>
#include <sstream>
#include <vector>

/** Description of a node in XML tree. Root node is the only node 
 * with NULL parent.*/
struct XMLNode {
   XMLNode* parent;
   std::multimap<std::string,XMLNode*> children;
   std::map<std::string,std::string> attributes;
   std::string value;
   
   XMLNode(XMLNode* parent);
   ~XMLNode();
};

class MuXML {
 public:
   MuXML();
   ~MuXML();

   template<typename T> bool addAttribute(XMLNode* node,const std::string& attribName,const T& attribValue);
   template<typename T> XMLNode* addNode(XMLNode* parent,const std::string& nodeName,const T& nodeValue);
   template<typename T> bool changeValue(XMLNode* node,const T& value);
   void clear();
   XMLNode* find(const std::string& nodeName,const XMLNode* node = NULL) const;
   XMLNode* find(const std::string& nodeName,const std::list<std::pair<std::string,std::string> >& attribs,const XMLNode* node=NULL) const;

   std::string getAttributeValue(const XMLNode* node,const std::string& attribName) const;
   void getAttributes(const XMLNode* node,std::map<std::string,std::string>& attribs) const;
   std::string getNodeValue(const XMLNode* node) const;
   XMLNode* getRoot() const;

   void print(std::ostream& out,const int& level=0,const XMLNode* node=NULL) const;

   bool read(std::istream& in,XMLNode* parent=NULL,const int& level=0,const char& currentChar=' ');
   
 private:
   
   XMLNode* root; /**< Pointer to root node.*/
};

template<typename T> bool MuXML::addAttribute(XMLNode* node,const std::string& attribName,const T& attribValue) {
   if (node == NULL) return false;
   // Add attribute through stringstream:
   std::stringstream ss;
   ss << attribValue;
   (node->attributes)[attribName] = ss.str();
   return true;
}

template<typename T> XMLNode* MuXML::addNode(XMLNode* parent,const std::string& nodeName,const T& nodeValue) {
   if (parent == NULL) return NULL;
   
   // Insert node:
   XMLNode* node = new XMLNode(parent);
   //(parent->children)[nodeName] = node;
   parent->children.insert(std::make_pair(nodeName,node));
   
   // Copy node value through stringstream:
   std::stringstream ss;
   ss << nodeValue;
   node->value = ss.str();   
   return node;
}

template<typename T> bool MuXML::changeValue(XMLNode* node,const T& value) {
   if (node == NULL) return false;
   
   // Change value through stringstream:
   std::stringstream ss;
   ss << value;
   node->value = ss.str();

   return true;
}

#endif

