/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#include <cstdlib>
#include <iostream>
#include "muxml.h"

using namespace std;

XMLNode::XMLNode(XMLNode* parent): parent(parent) { }

XMLNode::~XMLNode() {
   // Delete all children:
   for (multimap<string,XMLNode*>::iterator it=children.begin(); it!=children.end(); ++it) {
      delete it->second;
      it->second = NULL;
   }
}



MuXML::MuXML() {
   root = new XMLNode(NULL);
}

MuXML::~MuXML() {
   delete root;
   root = NULL;
}

void MuXML::clear() {
   delete root;
   root = new XMLNode(NULL);
}

XMLNode* MuXML::find(const std::string& nodeName,const XMLNode* node) const {
   if (node == NULL) node = root;
   // Recursive depth-first find. First check if child's name matches the searched name. If it does, return 
   // pointer to child. Otherwise search through child's children to see if it has the searched node.
   for (multimap<string,XMLNode*>::const_iterator it=node->children.begin(); it!=node->children.end(); ++it) {
      if (it->first == nodeName) return it->second;
      XMLNode* tmp = find(nodeName,it->second);
      if (tmp != NULL) return tmp;
   }
   return NULL;
}

XMLNode* MuXML::find(const std::string& nodeName,const std::list<std::pair<std::string,std::string> >& attribs,const XMLNode* node) const {
   if (node == NULL) node = root;
   list<pair<string,string> >::const_iterator jt;
   for (multimap<string,XMLNode*>::const_iterator it=node->children.begin(); it!=node->children.end(); ++it) {
      bool matchFound = true;
      
      if (it->first != nodeName) {
	 matchFound = false;
      }
      // Tag name matches, check that attributes match:
      if (matchFound == true) {
	 for (jt = attribs.begin(); jt!=attribs.end(); ++jt) {
	    map<string,string>::const_iterator tmp = it->second->attributes.find((*jt).first);
	    if (tmp == it->second->attributes.end()) {matchFound = false; break;} // attribute name was not found
	    if (tmp->second != (*jt).second) {matchFound = false; break;} // attribute value did not match
	 }
	 if (matchFound == true) {
	    return it->second;
	 }
      }
      // Recursively check children's nodes:
      XMLNode* tmp = find(nodeName,attribs,it->second);
      if (tmp != NULL) return tmp;
   }
   return NULL;
}

string MuXML::getAttributeValue(const XMLNode* node,const std::string& attribName) const {
   map<string,string>::const_iterator it=node->attributes.find(attribName);
   if (it == node->attributes.end()) return "";
   return it->second;
}

void MuXML::getAttributes(const XMLNode* node,std::map<std::string,std::string>& attribs) const {
   attribs = node->attributes;
}

string MuXML::getNodeValue(const XMLNode* node) const {
   return node->value;
}

XMLNode* MuXML::getRoot() const {return root;}

void MuXML::print(std::ostream& out,const int& level,const XMLNode* node) const {
   const int tab = 3;
   if (node == NULL) {
      node = root;
      //out << "XML TREE CONTENTS:" << endl;
   }
   for (multimap<string,XMLNode*>::const_iterator it=node->children.begin(); it!=node->children.end(); ++it) {
      // Indent
      for (int i=0; i<level; ++i) out << ' ';
      
      // Write child's name and its attributes:
      out << '<' << it->first;
      for (map<string,string>::const_iterator jt=it->second->attributes.begin(); jt!=it->second->attributes.end(); ++jt) {
	 out << ' ' << jt->first << "=\"" << jt->second << "\"";
      }
      
      // Write child's value:
      out << ">" << it->second->value;
      
      // Call print for the child:
      if (it->second->children.size() > 0) {
	 out << endl;
	 print(out,level+tab,it->second);
	 for (int i=0; i<level; ++i) out << ' ';
      }
      
      // Write child's end tag:
      out << "</" << it->first << '>' << endl;
   }
}

bool MuXML::read(std::istream& in,XMLNode* parent,const int& level,const char& currentChar) {
   in >> noskipws;
   if (parent == NULL) parent = root;
   
   int index = 0;
   bool success = true;
   char c = currentChar;
   char buffer[1024];
   do {
      // Skip empty chars:
      while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
      
      // Start to read tag name and its attributes:
      if (c == '<') {
	 in >> c;
	 while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	 
	 // Read end tag and return:
	 if (c == '/') {
	    while (in.eof() == false && c != '>') in >> c;
	    in >> c;
	    if (in.good() == false) success = false;
	    return success;
	 }
	 
	 index = 0;
	 while (c != ' ' && c != '>') {
	    buffer[index] = c;
	    ++index;
	    in >> c;
	    if (in.good() == false) {success = false; break;}
	 }
	 buffer[index] = '\0';
	 XMLNode* node = addNode(parent,buffer,"");
	 
	 // Remove empty spaces
	 while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	 if (in.eof() == true) {success = false; break;}
	 
	 // If next char is '>' the tag ends. Otherwise read attribute values:
	 if (c != '>') {
	    while (c != '>') {
	       while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	       
	       index = 0;
	       while (c != '=' && in.good() == true) {
		  buffer[index] = c;
		  ++index;
		  in >> c;
	       }
	       buffer[index] = '\0';
	       string attribName = buffer;
	       in >> c; // Forward from '='
	       in >> c; // Forward from '"'
	       index = 0;
	       while (c != '"' && in.good() == true) {
		  buffer[index] = c;
		  ++index;
		  in >> c;
	       }
	       in >> c;
	       if (in.good() == false) {success = false; break;}
	       buffer[index] = '\0';
	       string attribValue = buffer;
	       addAttribute(node,attribName,attribValue);
	    }
	    in >> c;
	 } else {
	    in >> c;
	 }
	 while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	 
	 // Read tag's value:
	 index = 0;
	 while ((c != ' ' && c != '\t' && c != '\n' && c != '<') && in.good() == true) {
	    buffer[index] = c;
	    ++index;
	    in >> c;
	 }
	 if (in.good() == false) {success = false; break;}
	 while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	 
	 buffer[index] = '\0';
	 changeValue(node,buffer);

	 if (c == '<') {
	    read(in,node,level+1,c);
	    in >> c;
	 }
	 while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
      }
      
      if (in.good() == false) return false;
   } while (success == true);
   return true;
}







