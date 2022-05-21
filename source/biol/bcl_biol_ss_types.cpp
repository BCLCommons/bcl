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
#include "biol/bcl_biol_ss_types.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_angle.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief constructor all SSTypes
    SSTypes::SSTypes() :
//                                                               IsStructured                  Radial Extent                   AnglePerTurn                  RiseInZPerRes                            Phi                            Psi                     FragLength            ContactWindowRadius            ThreeStatePrediction      backbone phi range      backbone psi range
      HELIX(              AddEnum( "HELIX"           , SSTypeData( 'H',  true,                         4.240,  math::Angle::Radian( -100.0),                       1.50247,                     -0.996767,                      -0.80718,                             5,                             4, linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>( math::Angle::Radian( -135.0), math::Angle::Radian( -25.0)), math::Range< double>( math::Angle::Radian( -70.0), math::Angle::Radian(  20.0))))),
      STRAND(             AddEnum( "STRAND"          , SSTypeData( 'E',  true,                         3.275,  math::Angle::Radian( -180.0),                         3.425,                      -2.39835,                       2.35502,                             3,                             2, linal::Vector3D( 0.0, 1.0, 0.0), math::Range< double>( math::Angle::Radian( -180.0), math::Angle::Radian( -35.0)), math::Range< double>( math::Angle::Radian(  25.0), math::Angle::Radian( 180.0))))),
      COIL(               AddEnum( "COIL"            , SSTypeData( 'C', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 0.0, 0.0, 1.0), math::Range< double>(), math::Range< double>()))),
      e_HelixRightOmega(  AddEnum( "HelixRightOmega" , SSTypeData( '2', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed omega
      e_HelixRightPi(     AddEnum( "HelixRightPi"    , SSTypeData( '3', false, util::GetUndefined< double>(),  math::Angle::Radian(  -87.0),                          1.15, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed pi
      e_HelixRightGamma(  AddEnum( "HelixRightGamma" , SSTypeData( '4', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed gamma
      e_HelixRight310(    AddEnum( "HelixRight310"   , SSTypeData( '5', false, util::GetUndefined< double>(),  math::Angle::Radian( -120.0),                           2.0, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix right handed 3-10
      e_HelixLeftAlpha(   AddEnum( "HelixLeftAlpha"  , SSTypeData( '6', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>( math::Angle::Radian(   45.0), math::Angle::Radian(  65.0)), math::Range< double>( math::Angle::Radian(  15.0), math::Angle::Radian( 100.0))))), // helix left handed alpha
      e_HelixLeftOmega(   AddEnum( "HelixLeftOmega"  , SSTypeData( '7', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix left handed omega
      e_HelixLeftGamma(   AddEnum( "HelixLeftGamma"  , SSTypeData( '8', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix left handed gamma
      e_Helix27Ribbon(    AddEnum( "Helix27Ribbon"   , SSTypeData( '9', false, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>(), math::Range< double>()))), // helix 2 - 7 ribbon
      e_HelixPolyProline( AddEnum( "HelixPolyProline", SSTypeData( '0', false, util::GetUndefined< double>(),  math::Angle::Radian(  120.0),                           3.1, util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), util::GetUndefined< size_t>(), linal::Vector3D( 1.0, 0.0, 0.0), math::Range< double>( math::Angle::Radian(  -80.0), math::Angle::Radian( -70.0)), math::Range< double>( math::Angle::Radian( 145.0), math::Angle::Radian( 155.0)))))  // helix poly proline
    {
      m_HelixTypes = storage::Set< SSType>( Begin() + e_HelixRightOmega, Begin() + e_HelixPolyProline + 1);
      m_HelixTypes.Insert( HELIX);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the reduced (3-state) types
    const storage::Vector< SSType> &SSTypes::GetReducedTypes() const
    {
      static storage::Vector< SSType> s_reduced( storage::Vector< SSType>::Create( HELIX, STRAND, COIL));
      return s_reduced;
    }

    //! @brief Get secondary structure type from provided ONE_LETTER_CODE
    //! @param ONE_LETTER_CODE one letter code for SSType
    //! @return secondary structure type from provided ONE_LETTER_CODE
    const SSType &SSTypes::SSTypeFromOneLetterCode( const char ONE_LETTER_CODE) const
    {
      // iterate over all ss types
      for
      (
        const_iterator type_itr( Begin()), type_itr_end( End());
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // if given ONE_LETTER_CODE matches this one letter code
        if( ( *type_itr)->GetOneLetterCode() == ONE_LETTER_CODE)
        {
          // return
          return *type_itr;
        }
      }

      // if no match is found return undefined
      return GetSSTypes().e_Undefined;
    }

    //! @brief secondary structure from phi and psi angle
    //! @param PHI phi backbone angle
    //! @param PSI psi backbone angle
    //! @return most likely SSType - if non was found, COIL
    const SSType &SSTypes::SSTypeFromPhiPsi( const double PHI, const double PSI) const
    {
      for( const_iterator sstype_itr( Begin()), sstype_itr_end( End()); sstype_itr != sstype_itr_end; ++sstype_itr)
      {
        if( ( *sstype_itr)->GetBackbonePhiRange().IsWithin( PHI) && ( *sstype_itr)->GetBackbonePsiRange().IsWithin( PSI))
        {
          return *sstype_itr;
        }
      }

      return COIL;
    }

    //! @brief converts given pdb helix class into SSType
    //! @see http://www.wwpdb.org/documentation/format32/sect5.html#HELIX
    //! @param HELIX_CLASS as found in HelixClass entry in pdb HELIX line 1-10
    //! @return helix SSType for given class, e_Undefined if HELIX_CLASS is not in [1,10] range
    const SSType &SSTypes::SSTypeFromPDBHelixClass( const size_t HELIX_CLASS) const
    {
      // from: http://www.wwpdb.org/documentation/format32/sect5.html#HELIX
      // TYPE OF  HELIX                     (COLUMNS 39 - 40)
      // --------------------------------------------------------------
      // Right-handed alpha (default)                1
      // Right-handed omega                          2
      // Right-handed pi                             3
      // Right-handed gamma                          4
      // Right-handed 3 - 10                         5
      // Left-handed alpha                           6
      // Left-handed omega                           7
      // Left-handed gamma                           8
      // 2 - 7 ribbon/helix                          9
      // Polyproline                                10
      static const math::Range< size_t> s_valid_classes( 1, 10);

      if( !s_valid_classes.IsWithin( HELIX_CLASS))
      {
        return HELIX;
//        return e_Undefined;
      }

      // right handed alpha helix
      if( HELIX_CLASS == 1)
      {
        return HELIX;
      }

      // remaining helices
      return *( e_HelixRightOmega.GetIterator() + ( HELIX_CLASS - s_AlphaOmegaHelixOffset));
    }

    //! @brief converts a SSType to a pdb helix class
    //! @param SSTYPE the sstype to convert to helix class
    //! @return pdb helix class [1-10]; if SSTYPE is not a valid helix class, undefined
    size_t SSTypes::PDBHelixClassFromSSType( const SSType &SSTYPE) const
    {
      // alpha helix is 1
      if( SSTYPE == HELIX)
      {
        return 1;
      }

      // strand coil undefined are undefined
      if( !SSTYPE.IsDefined() || SSTYPE <= COIL)
      {
        return util::GetUndefined< size_t>();
      }

      // rest is directly related to ss type
      const size_t helix_class( SSTYPE.GetIndex() + 1 - s_AlphaOmegaHelixOffset);
      if( helix_class > 10)
      {
        return util::GetUndefined< size_t>();
      }

      return helix_class;
    }

    //! @brief set of all helix types
    //! @return Set containing all helix types
    const storage::Set< SSType> &SSTypes::GetHelixTypes() const
    {
      return m_HelixTypes;
    }

    //! @brief test if two sses are of similar type (alpha helix to 3-10 helix are considered similar, but strand would not be)
    //! @param SS_TYPE_LHS left hand side ss type
    //! @param SS_TYPE_RHS right hand side ss type
    //! @return true if the two types are considered similar
    bool SSTypes::AreSimilar( const SSType &SS_TYPE_LHS, const SSType &SS_TYPE_RHS) const
    {
      // identical type is similar
      if( SS_TYPE_LHS == SS_TYPE_RHS)
      {
        return true;
      }

      // consider right helices as similar
      if( SS_TYPE_LHS == HELIX)
      {
        if( SS_TYPE_RHS >= e_HelixRightOmega && SS_TYPE_RHS <= e_HelixRight310)
        {
          return true;
        }
      }
      if( SS_TYPE_RHS == HELIX)
      {
        if( SS_TYPE_LHS >= e_HelixRightOmega && SS_TYPE_LHS <= e_HelixRight310)
        {
          return true;
        }
      }

      return false;
    }

    //! @brief on access function for all SSTypes
    //! @return SSTypes enums
    const SSTypes &GetSSTypes()
    {
      return SSTypes::GetEnums();
    }

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< biol::SSTypeData, biol::SSTypes>;

  } // namespace util
} // namespace bcl
