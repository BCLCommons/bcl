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

#ifndef BCL_SSPRED_DSSP_STRIDE_H_
#define BCL_SSPRED_DSSP_STRIDE_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_sspred_dssp.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DsspStride
    //! @brief stores prediction for DsspStride, which merges non-coil predictions from these two methods
    //! The idea is that Stride does a great job of identifying the ends of SSEs, generally choosing more residues
    //! than DSSP, however, the problem is that Stride sometimes misses entire helices that are detected by DSSP.
    //! In benchmarks on both soluble and membrane proteins, about 2-3% of all Stride-declared Coiled residues are
    //! correctly identified as helix by DSSP
    //! simple script that runs and merges these two methods predictions according to the following matrix:
    //!               | DSSP Helix | DSSP Strand | DSSP Coil |
    //! Stride Helix  | Helix      | Coil        | Helix     |
    //! Stride Strand | Coil       | Strand      | Strand    |
    //! Stride Coil   | Helix      | Strand      | Coil      |
    //!
    //! @see @link example_sspred_dssp_stride.cpp @endlink
    //! @author mendenjl
    //! @date Aug 19, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DsspStride :
      public Dssp
    {

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DsspStride();

      //! @brief constructor from an sse type
      //! @param SS_TYPE type of sse
      DsspStride( const biol::SSType &SS_TYPE);

      //! @brief Clone function
      //! @return pointer to new DsspStride
      DsspStride *Clone() const;

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

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read secondary structure predictions for given sequence from the provided ISTREAM
      //! @param ISTREAM input stream
      //! @param AA_SEQUENCE AASequence into which sspredictions will be read
      //! @return std::istream which was read from
      std::istream &ReadPredictionsForAASequence( std::istream &ISTREAM, biol::AASequence &AA_SEQUENCE) const;

    }; // class DsspStride

  } // namespace sspred
} // namespace bcl

#endif // BCL_SSPRED_DSSP_STRIDE_H_
