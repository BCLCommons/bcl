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

#ifndef BCL_SSPRED_PALSSE_H_
#define BCL_SSPRED_PALSSE_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_method_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Palsse
    //! @brief stores prediction for Palsse
    //! @details stores prediction for Palsse, which is a
    //!
    //! @see @link example_sspred_palsse.cpp @endlink
    //! @author mendenjl
    //! @date Aug 19, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Palsse :
      public MethodInterface
    {

    private:

    //////////
    // data //
    //////////

      //! native SSType defined by the Palsse algorithm
      biol::SSType m_SSType;

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Palsse();

      //! @brief constructor from an sse type
      //! @param SS_TYPE type of sse
      Palsse( const biol::SSType &SS_TYPE);

      //! @brief Clone function
      //! @return pointer to new Palsse
      Palsse *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get file extension associated with this Method
      //! @return file extension associated with this Method
      const std::string &GetFileExtension() const;

      //! @brief get whether this method determined the secondary structure / membrane environment from the structure
      //! @return true if this method determined the secondary structure / membrane environment from the structure
      bool GetIsDeterminedFromSturcture() const
      {
        return true;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get three state, environment independent secondary structure prediction
      //! @return three state, environment independent secondary structure prediction
      linal::Vector3D GetThreeStatePrediction() const;

      //! @brief get nine state secondary structure prediction ( 3 SSTypes for all 3 EnvironmentTypes)
      //! @return three state secondary structure prediction
      linal::Matrix< double> GetNineStatePrediction() const;

      //! @brief read secondary structure predictions for given amino acid from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AMINO_ACID amino acid into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAA( std::istream &ISTREAM, biol::AABase &AMINO_ACID) const;

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const;

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

    }; // class Palsse

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_PALSSE_H_
