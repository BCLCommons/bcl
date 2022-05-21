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

#ifndef BCL_RESTRAINT_HANDLER_DATA_SET_PAIRWISE_INTERFACE_H_
#define BCL_RESTRAINT_HANDLER_DATA_SET_PAIRWISE_INTERFACE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerDataSetPairwiseInterface
    //! @brief Provides functionality for reading and writing pairwise data sets in specific formats
    //! @details none
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date Oct 28, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HandlerDataSetPairwiseInterface :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reads data set from a stream
      //! @param ISTREAM the stream the data set will be read from
      //! @return stream that the data set was read from
      virtual std::istream &ReadDataSetPairwise( std::istream &ISTREAM) = 0;

      //! @brief provides the data set that the handler created
      //! @return const reference to an data set object
      virtual const DataSetPairwise &GetDataSetPairwise() const = 0;

      //! @brief writes restraints to a stream
      //! @param OSTREAM the stream the data set will be written to
      //! @param DATA_SET the data set that will be written to the stream
      //! @return stream the data set were written to
      virtual std::ostream &WriteDataSetPairwise( std::ostream &OSTREAM, const DataSetPairwise &DATA_SET) const = 0;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class HandlerDataSetPairwiseInterface

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_HANDLER_DATA_SET_PAIRWISE_INTERFACE_H_ 
