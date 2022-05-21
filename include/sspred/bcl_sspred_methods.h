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

#ifndef BCL_SSPRED_METHODS_H_
#define BCL_SSPRED_METHODS_H_

// include the namespace header
#include "bcl_sspred.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Methods
    //! @brief Enumerator class for different secondary structure methods
    //! @details This enumerator class different secondary structure methods stores a ShPtr to MethodInterface for each enum.
    //!
    //! @see @link example_sspred_methods.cpp @endlink
    //! @author karakam
    //! @date Jun 3, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Methods :
      public util::Enumerate< util::ShPtr< MethodInterface>, Methods>
    {
      friend class util::Enumerate< util::ShPtr< MethodInterface>, Methods>;

    public:

    //////////
    // data //
    //////////

      Method e_PDB;             //!< PDB
      Method e_PSIPRED;         //!< PSIPRED
      Method e_JUFO;            //!< JUFO
      Method e_JUFO9D;          //!< JUFO9D
      Method e_SAM;             //!< SAM
      Method e_PROFphd;         //!< PROFphd
      Method e_TMHMM;           //!< TMHMM
      Method e_TMMOD;           //!< TMMOD
      Method e_B2TMPRED;        //!< B2TMPRED
      Method e_PROFTMB;         //!< PROFTMB
      Method e_CONPRED;         //!< CONPRED
      Method e_TALOS;           //!< TALOS
      Method e_OCTOPUS;         //!< OCTOPUS
      Method e_BOCTOPUS;        //!< BOCTOPUS
      Method e_TMBETANET;       //!< TMBETANET
      Method e_PARTIFOLD;       //!< PARTIFOLD
      Method e_MASP;            //!< Membrane association & Secondary structure prediction
      Method e_Stride;          //!< Stride analysis method (http://webclu.bio.wzw.tum.de/stride/)
      Method e_DSSP;            //!< DSSP analysis method (http://swift.cmbi.ru.nl/gv/dssp/)
      Method e_StrideDSSP;      //!< Merges predictions of Stride with those of DSSP, leveraging the strengths of both
      Method e_PALSSE;          //!< PALSSE vector based analysis method
      Method e_MAHSSMI;         //!< Membrane Aware Hybrid Secondary Structure & Membrane topology Identification
      Method e_CIPhiPsi;        //!< Context-insensitive phi-psi based secondary structure analysis
      Method e_Kaksi;           //!< KAKSI analysis method

      //! @brief return command line flag for specifying which ss prediction methods to read in
      //! @return command line flag for specifying which ss prediction methods to read in
      static util::ShPtr< command::FlagInterface> &GetFlagReadSSPredictions();

      //! @brief get vector of SS methods
      //! @return vector of sspred::Methods which were given over the command line
      static storage::Set< Method> GetCommandLineMethods();

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor that constructs all Methods
      Methods();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

    }; // class Methods

    //! @brief construct on access function for all Methods
    //! @return reference to only instances of Methods
    BCL_API const Methods &GetMethods();

  } // namespace sspred

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< sspred::MethodInterface>, sspred::Methods>;

  } // namespace util
} // namespace bcl

#endif // BCL_SSPRED_METHODS_H_
