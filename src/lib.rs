mod tests;
mod spooles_binding;
use crate::spooles_binding::
{
    SPOOLES_REAL, SPOOLES_SYMMETRIC, SPOOLES_NONSYMMETRIC, INPMTX_BY_ROWS, INPMTX_BY_VECTORS, 
    INPMTX_BY_CHEVRONS, NO_LOCK, FRONTMTX_DENSE_FRONTS, SPOOLES_NO_PIVOTING, SPOOLES_PIVOTING, IV, IVL, FILE,
    InpMtx_new, InpMtx_init, InpMtx_inputRealEntry, InpMtx_changeStorageMode, DenseMtx_new, DenseMtx_init, 
    DenseMtx_setRealEntry, Graph_new, InpMtx_fullAdjacency, IVL_tsize, Graph_init2, orderViaMMD, ETree_oldToNewVtxPerm, 
    IV_entries, ETree_newToOldVtxPerm, ETree_permuteVertices, InpMtx_permute, InpMtx_mapToUpperTriangle, 
    InpMtx_changeCoordType, DenseMtx_permuteRows, SymbFac_initFromInpMtx, FrontMtx_new, SubMtxManager_new,
    SubMtxManager_init, FrontMtx_init, ChvManager_new, ChvManager_init, DVfill, IVfill, FrontMtx_factorInpMtx,
    FrontMtx_postProcess, DenseMtx_zero, FrontMtx_solve, DenseMtx_entries, FrontMtx_free, DenseMtx_free,
    IV_free, InpMtx_free, ETree_free, IVL_free, SubMtxManager_free, Graph_free,
};


pub enum SymmetryFlag
{
    Symmetric,
    NonSymmetric,
}


pub enum PivotingFlag
{
    Pivoting,
    NoPivoting,
}


/// Serial solution of AX = B using LU factorization.
pub fn solve_using_lu(n_row: i32, a: Vec<(i32, i32, f64)>, b: Vec<(i32, f64)>, symmetry_flag: SymmetryFlag, 
    pivoting_flag: PivotingFlag) -> Result<Vec<f64>, String>
{
    let symmetry_flag = match symmetry_flag 
        {
            SymmetryFlag::Symmetric => SPOOLES_SYMMETRIC,
            SymmetryFlag::NonSymmetric => SPOOLES_NONSYMMETRIC,
        };

    let pivoting_flag = match pivoting_flag
        {
            PivotingFlag::Pivoting => SPOOLES_PIVOTING,
            PivotingFlag::NoPivoting => SPOOLES_NO_PIVOTING,
        };
    
    let n_ent = a.len() as i32;

    let m_type = SPOOLES_REAL;
    let n_eqns = n_row;

    let n_rhs = 1;
    let seed = 1;
    let msg_lvl = 0;
    let msg_file: *mut FILE = std::ptr::null_mut();
    let mut cpus = [0f64; 10];
    let mut stats = [0; 20];
    let droptol = 0.0; 
    let tau = 100.0;
    let mut error = 0;

    let mtx_a = unsafe{ InpMtx_new() };
    unsafe 
    { 
        InpMtx_init(mtx_a, INPMTX_BY_ROWS as i32, m_type as i32, n_ent, 
            n_eqns); 
    }
    for (i_row, j_col, a_value) in a.iter()
    {
        unsafe { InpMtx_inputRealEntry(mtx_a, *i_row, *j_col, *a_value); }
    }
    unsafe { InpMtx_changeStorageMode(mtx_a, INPMTX_BY_VECTORS as i32); }

    let mtx_y = unsafe { DenseMtx_new() };
    unsafe 
    { 
        DenseMtx_init(mtx_y, m_type as i32, 0, 0, n_eqns, n_rhs, 1, n_eqns); 
    }
    for (i_row, y_value) in b.iter()
    {
        unsafe { DenseMtx_setRealEntry(mtx_y, *i_row, 0, *y_value); }
    }

    let graph = unsafe { Graph_new() };
    let adj_ivl = unsafe { InpMtx_fullAdjacency(mtx_a) };
    let n_edges = unsafe { IVL_tsize(adj_ivl) };
    let v_vghts: *mut i32 = std::ptr::null_mut();
    let e_wght_ivl: *mut IVL = std::ptr::null_mut();
    unsafe 
    {
        Graph_init2(graph, 0, n_eqns, 0, n_edges, n_eqns, n_edges, 
            adj_ivl, v_vghts, e_wght_ivl);
    }

    let front_e_tree = unsafe { orderViaMMD(graph, seed, msg_lvl, msg_file) };
    let old_to_new_iv = unsafe { ETree_oldToNewVtxPerm(front_e_tree) };
    let old_to_new = unsafe { IV_entries(old_to_new_iv) };
    let new_to_old_iv = unsafe { ETree_newToOldVtxPerm(front_e_tree) };
    unsafe 
    {
        ETree_permuteVertices(front_e_tree, new_to_old_iv);
        InpMtx_permute(mtx_a, old_to_new, old_to_new);
        if symmetry_flag == SPOOLES_SYMMETRIC
        {
            InpMtx_mapToUpperTriangle(mtx_a);
        }
        InpMtx_changeCoordType(mtx_a, INPMTX_BY_CHEVRONS as i32);
        InpMtx_changeStorageMode(mtx_a, INPMTX_BY_VECTORS as i32);
        DenseMtx_permuteRows(mtx_y, old_to_new_iv);
    }
    let symb_fac_ivl = unsafe { SymbFac_initFromInpMtx(front_e_tree, mtx_a) };

    let front_mtx = unsafe { FrontMtx_new() };
    let mtx_manager = unsafe { SubMtxManager_new() };
    let owners_iv: *mut IV = std::ptr::null_mut();
    unsafe
    {
        SubMtxManager_init(mtx_manager, NO_LOCK as i32, 0);
        FrontMtx_init(front_mtx, front_e_tree, symb_fac_ivl, m_type as i32, 
            symmetry_flag as i32, FRONTMTX_DENSE_FRONTS as i32, 
            pivoting_flag as i32, NO_LOCK as i32, 0, owners_iv, mtx_manager, 
            msg_lvl, msg_file);
    }

    let chv_manager = unsafe { ChvManager_new() };
    unsafe 
    {
        ChvManager_init(chv_manager, NO_LOCK as i32, 1);
        DVfill(10, cpus.as_mut_ptr(), 0.0);
        IVfill(20, stats.as_mut_ptr(), 0);
    }
    let root_chv = unsafe { FrontMtx_factorInpMtx(front_mtx, mtx_a, tau, droptol, 
        chv_manager, &mut error, cpus.as_mut_ptr(), stats.as_mut_ptr(), msg_lvl,
        msg_file) };

    if root_chv != std::ptr::null_mut()
    {
        return Err(String::from("Matrix found to be singular!"));
    }

    if error >= 0
    {
        return Err(format!("Error encountered at front {error}!"));
    }

    unsafe { FrontMtx_postProcess(front_mtx, msg_lvl, msg_file); }

    let mtx_x = unsafe { DenseMtx_new() };
    unsafe
    {
        DenseMtx_init(mtx_x, m_type as i32, 0, 0, n_eqns, n_rhs, 1, n_eqns);
        DenseMtx_zero(mtx_x);
        FrontMtx_solve(front_mtx, mtx_x, mtx_y, mtx_manager, cpus.as_mut_ptr(), 
            msg_lvl, msg_file);
    }

    unsafe { DenseMtx_permuteRows(mtx_x, new_to_old_iv); }

    let result = unsafe { std::slice::from_raw_parts(DenseMtx_entries(mtx_x), n_row as usize) }
        .to_vec();

    unsafe 
    {
        FrontMtx_free(front_mtx);
        DenseMtx_free(mtx_x);
        DenseMtx_free(mtx_y);
        IV_free(new_to_old_iv);
        IV_free(old_to_new_iv);
        InpMtx_free(mtx_a);
        ETree_free(front_e_tree);
        IVL_free(symb_fac_ivl);
        SubMtxManager_free(mtx_manager);
        Graph_free(graph);
    }

    Ok(result)
}
