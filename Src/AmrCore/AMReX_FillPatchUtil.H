#ifndef AMREX_FillPatchUtil_H_
#define AMREX_FillPatchUtil_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MFInterpolater.H>
#include <AMReX_Array.H>
#include <AMReX_Utility.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBInterpolater.H>
#include <AMReX_EBMFInterpolater.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cmath>
#include <limits>

namespace amrex
{

    template <typename FAB>
    struct NullInterpHook
    {
        void operator() (FAB& /*fab*/, const Box& /*bx*/, int /*icomp*/, int /*ncomp*/) const {}

        void operator() (Array<FAB*, AMREX_SPACEDIM> /*fab*/, const Box& /*bx*/, int /*icomp*/, int /*ncomp*/) const {}
    };

    template <typename Interp>
    bool ProperlyNested (const IntVect& ratio, const IntVect& blocking_factor, int ngrow,
                         const IndexType& boxType, Interp* mapper);

    
    //SWAPNIL: This is the main function for FillPatchSingleLevel. 
    //That's because this one is ultimately called by the next function later.
    /*****************   MAIN **************************************************/
    template <typename MF, typename BC>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (MF& mf, IntVect const& nghost, Real time,
                          const Vector<MF*>& smf, const Vector<Real>& stime,
                          int scomp, int dcomp, int ncomp,
                          const Geometry& geom,
                          BC& physbcf, int bcfcomp);  
    /****************************************************************************/
    
    /*****************   DEPENDENT (USED IN CARPETX) **************************************************/
    //No need to specify the vector of ghost-zones in this one, taken care of from mf ultimately
    template <typename MF, typename BC>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (MF& mf, Real time,
                          const Vector<MF*>& smf, const Vector<Real>& stime,
                          int scomp, int dcomp, int ncomp,
                          const Geometry& geom,
                          BC& physbcf, int bcfcomp); //SWAPNIL: This is used in AmrCoreAdv.cpp
    /****************************************************************************/
     
    //------------------------- My version: take vector of MFs as input ----------------------------
    template <typename MF, typename BC>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (amrex::Vector<MF*> const& vec_mf, amrex::Vector<amrex::IntVect>& vec_nghost, 
                          const amrex::Vector<Real>& vec_time,
                          const amrex::Vector<amrex::Vector<MF*>>& vec_vec_smf, 
                          const Vector<Vector<Real>>& vec_vec_stime, amrex::Vector<int>& vec_scomp, 
                          const amrex::Vector<int>& vec_dcomp, const amrex::Vector<int>& vec_ncomp,
                          const amrex::Vector<Geometry>& vec_geom,
                          amrex::Vector<BC>& vec_physbcf, const amrex::Vector<int>& vec_bcfcomp);
    
    template <typename MF, typename BC>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchSingleLevel (amrex::Vector<MF*>const& vec_mf, const amrex::Vector<Real>& vec_time,
                          const amrex::Vector<amrex::Vector<MF*>>& vec_vec_smf, 
                          const Vector<Vector<Real>>& vec_vec_stime,
                          const amrex::Vector<int>& vec_scomp, const amrex::Vector<int>& vec_dcomp, 
                          const amrex::Vector<int>& vec_ncomp,
                          const amrex::Vector<Geometry>& vec_geom,
                          amrex::Vector<BC>& vec_physbcf, const amrex::Vector<int>& vec_bcfcomp);
    //------------------------------------------------------------------------------------------------
    
    
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, IntVect const& nghost, Real time,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    /*****************   (USED IN CARPETX) **************************************************/
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, Real time,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});
    //****************************************************************************
    
    //---------------------   MY VERSION WITH VEC -----------------------------------------------
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (amrex::Vector<MF*>& vec_mf, amrex::Vector<Real>& vec_time,
                        const amrex::Vector<amrex::Vector<MF*>>& vec_vec_cmf, 
                        const Vector<Vector<Real>>& vec_vec_ct,
                        const amrex::Vector<amrex::Vector<MF*>>& vec_vec_fmf,
                        const Vector<Vector<Real>>& vec_vec_ft,
                        amrex::Vector<int>& vec_scomp, amrex::Vector<int>& vec_dcomp, 
                        amrex::Vector<int>& vec_ncomp,
                        const amrex::Vector<Geometry>& vec_cgeom,
                        const amrex::Vector<Geometry>& vec_fgeom,
                        amrex::Vector<BC>& vec_cbc, amrex::Vector<int>& vec_cbccomp,
                        amrex::Vector<BC>& vec_fbc, amrex::Vector<int>& vec_fbccomp,
                        const amrex::Vector<IntVect>& vec_ratio,
                        amrex::Vector<Interp*>& vec_mapper,                           
                        const amrex::Vector<amrex::Vector<BCRec>>& vec_vec_bcs,
                        amrex::Vector<int>& vec_bcscomp,
                        const amrex::Vector<PreInterpHook>& vec_pre_interp = {{}},
                        const amrex::Vector<PostInterpHook>& vec_post_interp = {{}});
    //------------------------------------------------------------------------------------------------
    
    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        Array<BC, AMREX_SPACEDIM>& cbc, const Array<int, AMREX_SPACEDIM>& cbccomp,
                        Array<BC, AMREX_SPACEDIM>& fbc, const Array<int, AMREX_SPACEDIM>& fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, const Array<int, AMREX_SPACEDIM>& bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                        Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (Array<MF*, AMREX_SPACEDIM> const& mf, Real time,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& cmf, const Vector<Real>& ct,
                        const Vector<Array<MF*, AMREX_SPACEDIM> >& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                        Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp = {},
                        const PostInterpHook& post_interp = {});

#ifdef AMREX_USE_EB
    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, IntVect const& nghost, Real time,
                        const EB2::IndexSpace& index_space,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp,
                        const PostInterpHook& post_interp);

    template <typename MF, typename BC, typename Interp, typename PreInterpHook, typename PostInterpHook>
    std::enable_if_t<IsFabArray<MF>::value>
    FillPatchTwoLevels (MF& mf, Real time,
                        const EB2::IndexSpace& index_space,
                        const Vector<MF*>& cmf, const Vector<Real>& ct,
                        const Vector<MF*>& fmf, const Vector<Real>& ft,
                        int scomp, int dcomp, int ncomp,
                        const Geometry& cgeom, const Geometry& fgeom,
                        BC& cbc, int cbccomp,
                        BC& fbc, int fbccomp,
                        const IntVect& ratio,
                        Interp* mapper,
                        const Vector<BCRec>& bcs, int bcscomp,
                        const PreInterpHook& pre_interp,
                        const PostInterpHook& post_interp);
#endif

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (MF& mf, Real time,
                           const MF& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           BC& cbc, int cbccomp,
                           BC& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Vector<BCRec>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (MF& mf, IntVect const& nghost, Real time,
                           const MF& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           BC& cbc, int cbccomp,
                           BC& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Vector<BCRec>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM> const& mf, Real time,
                           const Array<MF*, AMREX_SPACEDIM>& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                           Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

    template <typename MF, typename BC, typename Interp,
              typename PreInterpHook=NullInterpHook<typename MF::FABType::value_type>,
              typename PostInterpHook=NullInterpHook<typename MF::FABType::value_type> >
    std::enable_if_t<IsFabArray<MF>::value>
    InterpFromCoarseLevel (Array<MF*, AMREX_SPACEDIM> const& mf, IntVect const& nghost, Real time,
                           const Array<MF*, AMREX_SPACEDIM>& cmf, int scomp, int dcomp, int ncomp,
                           const Geometry& cgeom, const Geometry& fgeom,
                           Array<BC, AMREX_SPACEDIM>& cbc, int cbccomp,
                           Array<BC, AMREX_SPACEDIM>& fbc, int fbccomp,
                           const IntVect& ratio,
                           Interp* mapper,
                           const Array<Vector<BCRec>, AMREX_SPACEDIM>& bcs, int bcscomp,
                           const PreInterpHook& pre_interp = {},
                           const PostInterpHook& post_interp = {});

#ifndef BL_NO_FORT
    enum InterpEM_t { InterpE, InterpB};

    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const Array<MultiFab,AMREX_SPACEDIM>& crse,
                                     Array<MultiFab,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio);

    void InterpCrseFineBndryEMfield (InterpEM_t interp_type,
                                     const Array<MultiFab const*,AMREX_SPACEDIM>& crse,
                                     const Array<MultiFab*,AMREX_SPACEDIM>& fine,
                                     const Geometry& cgeom, const Geometry& fgeom,
                                     int ref_ratio);
#endif
}

#include <AMReX_FillPatchUtil_I.H>

#endif
