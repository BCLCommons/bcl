// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_APP_SCATTERPLOT_FROM_TABLES_H_
#define BCL_APP_SCATTERPLOT_FROM_TABLES_H_

// include the namespace header

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScatterplotFromTables
    //! @brief This app generates scatter plot gnuplot files from scoring tables
    //! @details This app takes a scoring table output from app::FoldAnalysis and generates a gnuplot file for a
    //! scatter plot
    //!
    //! @see @link example_app_scatterplot_from_tables.cpp @endlink
    //! @author alexanns, weinerbe
    //! @date Mar 7, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScatterplotFromTables :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ScatterplotFromTables();

      //! @brief Clone function
      //! @return pointer to new ScatterplotFromTables
      ScatterplotFromTables *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the Command object
      util::ShPtr< command::Command> InitializeCommand() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType ScatterplotFromTables_Instance;

    }; // class ScatterplotFromTables

  } // namespace app
} // namespace bcl

#endif // BCL_APP_SCATTERPLOT_FROM_TABLES_H_
