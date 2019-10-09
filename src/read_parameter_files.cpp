/* ----------------------------------------------------------------------
    Copyright (C) 2019 Nicholas Huebner Julian <njulian@ucla.edu>
---------------------------------------------------------------------- */
// File: read_parameter_files.cpp

#ifndef READ_PARAMETER_FILES_CPP
#define READ_PARAMETER_FILES_CPP

#include "read_parameter_files.hpp"

int SPF_NS::read_state_file_list(
               const std::string& list_file_name,
               std::vector<std::string>& state_file_names
               )
{
   std::ifstream list_file( list_file_name.c_str() );
   if ( list_file.is_open() )
   {
      std::string state_file_name;
      std::string list_file_line;
      size_t first_pos; 
      // iterate over the lines in the file
      while (std::getline( list_file, list_file_line ) && list_file.good())
      {
         // skipping lines of whitespace
         first_pos = list_file_line.find_first_not_of(" \t");
         while ( first_pos == std::string::npos )
         {
            std::getline( list_file, list_file_line );
            first_pos = list_file_line.find_first_not_of(" \t");
         }
         std::istringstream list_file_line_stream( list_file_line );
         list_file_line_stream >> state_file_name;

         state_file_names.push_back( state_file_name );
      }
      //if ( ! list_file.good() )
      //{
      //   std::cerr << "Error, state of opened file "
      //      << list_file_name
      //      << ".good() is false" << std::endl;
      //   list_file.close();
      //   std::cout << "state_file_names : ";
      //   for ( std::vector<std::string>::const_iterator itr = state_file_names.begin(); itr != state_file_names.end(); ++itr)
      //      std::cout << std::endl << "   " << *itr;
      //   std::cout << std::endl;
      //   return EXIT_FAILURE;
      //}
   }
   else
   {
      std::cerr << "Error opening file " << list_file_name << std::endl;
      list_file.close();
      return EXIT_FAILURE;
   }
   list_file.close();
   return EXIT_SUCCESS;
}

//int SPF_NS::read_parameter_file(
//      )
//{
//   return EXIT_SUCCESS;
//}

#endif
