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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "sspred/bcl_sspred_methods.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "sspred/bcl_sspred_b2tmpred.h"
#include "sspred/bcl_sspred_boctopus.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "sspred/bcl_sspred_conpred.h"
#include "sspred/bcl_sspred_dssp_stride.h"
#include "sspred/bcl_sspred_jufo.h"
#include "sspred/bcl_sspred_jufo9d.h"
#include "sspred/bcl_sspred_kaksi.h"
#include "sspred/bcl_sspred_mahssmi.h"
#include "sspred/bcl_sspred_masp.h"
#include "sspred/bcl_sspred_octopus.h"
#include "sspred/bcl_sspred_palsse.h"
#include "sspred/bcl_sspred_partifold.h"
#include "sspred/bcl_sspred_pdb.h"
#include "sspred/bcl_sspred_profphd.h"
#include "sspred/bcl_sspred_proftmb.h"
#include "sspred/bcl_sspred_psipred.h"
#include "sspred/bcl_sspred_sam.h"
#include "sspred/bcl_sspred_stride.h"
#include "sspred/bcl_sspred_talos.h"
#include "sspred/bcl_sspred_tmbetanet.h"
#include "sspred/bcl_sspred_tmhmm.h"
#include "sspred/bcl_sspred_tmmod.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {
  //////////
  // data //
  //////////

    //! @brief return command line flag for specifying which ss prediction methods to read in
    //! @return command line flag for specifying which ss prediction methods to read in
    util::ShPtr< command::FlagInterface> &Methods::GetFlagReadSSPredictions()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "sspred",
          "\tFlag for reading in secondary structure predictions from specified methods",
          command::Parameter
          (
            "ssmethod",
            "\tmethodname for input secondary structure",
            command::ParameterCheckEnumerate< Methods>(),
            GetMethods().e_JUFO.GetName()
          ),
          0, GetMethods().GetEnumCount()
        )
      );
      // end
      return s_flag;
    }

    //! @brief get vector of SS methods
    //! @return vector of sspred::Methods which were given over the command line
    storage::Set< Method> Methods::GetCommandLineMethods()
    {
      return GetFlagReadSSPredictions()->GetObjectSet< Method>();
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor that constructs all Methods
    Methods::Methods() :
      e_PDB(           AddEnum( "PDB",           util::ShPtr< MethodInterface>( new PDB()))),
      e_PSIPRED(       AddEnum( "PSIPRED",       util::ShPtr< MethodInterface>( new PSIPRED()))),
      e_JUFO(          AddEnum( "JUFO",          util::ShPtr< MethodInterface>( new JUFO()))),
      e_JUFO9D(        AddEnum( "JUFO9D",        util::ShPtr< MethodInterface>( new JUFO9D()))),
      e_SAM(           AddEnum( "SAM",           util::ShPtr< MethodInterface>( new SAM()))),
      e_PROFphd(       AddEnum( "PROFphd",       util::ShPtr< MethodInterface>( new PROFphd()))),
      e_TMHMM(         AddEnum( "TMHMM",         util::ShPtr< MethodInterface>( new TMHMM()))),
      e_TMMOD(         AddEnum( "TMMOD",         util::ShPtr< MethodInterface>( new TMMOD()))),
      e_B2TMPRED(      AddEnum( "B2TMPRED",      util::ShPtr< MethodInterface>( new B2TMPRED()))),
      e_PROFTMB(       AddEnum( "PROFTMB",       util::ShPtr< MethodInterface>( new PROFTMB()))),
      e_CONPRED(       AddEnum( "CONPRED",       util::ShPtr< MethodInterface>( new CONPRED()))),
      e_TALOS(         AddEnum( "TALOS",         util::ShPtr< MethodInterface>( new TALOS()))),
      e_OCTOPUS(       AddEnum( "OCTOPUS",       util::ShPtr< MethodInterface>( new OCTOPUS()))),
      e_BOCTOPUS(      AddEnum( "BOCTOPUS",      util::ShPtr< MethodInterface>( new BOCTOPUS()))),
      e_TMBETANET(     AddEnum( "TMBETANET",     util::ShPtr< MethodInterface>( new TMBETANET()))),
      e_PARTIFOLD(     AddEnum( "PARTIFOLD",     util::ShPtr< MethodInterface>( new PARTIFOLD()))),
      e_MASP(          AddEnum( "MASP",          util::ShPtr< MethodInterface>( new MASP()))),
      e_Stride(        AddEnum( "Stride",        util::ShPtr< MethodInterface>( new Stride()))),
      e_DSSP(          AddEnum( "DSSP",          util::ShPtr< MethodInterface>( new Dssp()))),
      e_StrideDSSP(    AddEnum( "DSSPStride",    util::ShPtr< MethodInterface>( new DsspStride()))),
      e_PALSSE(        AddEnum( "Palsse",        util::ShPtr< MethodInterface>( new Palsse()))),
      e_MAHSSMI(       AddEnum( "MAHSSMI",       util::ShPtr< MethodInterface>( new Mahssmi()))),
      e_CIPhiPsi(      AddEnum( "CIPhiPsi",      util::ShPtr< MethodInterface>( new CIPhiPsi()))),
      e_Kaksi(         AddEnum( "Kaksi",         util::ShPtr< MethodInterface>( new Kaksi())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Methods::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct on access function for all Methods
    //! @return reference to only instances of Methods
    const Methods &GetMethods()
    {
      return Methods::GetEnums();
    }

  } // namespace sspred

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< sspred::MethodInterface>, sspred::Methods>;

  } // namespace util
} // namespace bcl
