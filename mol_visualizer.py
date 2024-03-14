def create_view(view, 
                mol_true='',
                mol_cond='',
                mol_pred='', 
                hide_true=False,
                hide_cond=False, 
                hide_pred=False,
                surface=None,
                opacity=0.5,
                show_lines=False,
                ):

    view.addModel(str(mol_cond), 'pdb')
    if not hide_cond:
        for i, at in enumerate(mol_cond):
            default = {'cartoon': {'color': '#666666', "opacity": 0.7}}
            if show_lines:
                default['line'] = {}
            
            final_style = {**at.get('style', default), **default}
            view.setStyle({'model': -1, 'serial': i+1}, final_style)
        if surface:
            view.addSurface(surface, {'opacity': opacity, 'color':'#222222'})

    if not hide_true:
        view.addModel(str(mol_true), 'sdf')
        for i, at in enumerate(mol_true):
            default = {'stick': {'colorscheme': 'greenCarbon'}}
            view.setStyle({'model': -1, 'serial': i}, at.get("style", default))


    if not hide_pred:
        view.addModel(str(mol_pred), 'sdf')
        for i, at in enumerate(mol_pred):
            default = {'stick': {'colorscheme': 'lightBlueCarbon'}}
            style = at.get("style", default)
            view.setStyle({'model': -1, 'serial': i}, style)

    view.zoomTo()

def draw_line(view, at0, at1, kwargs={}):
    view.addLine({'start': {'x': at0["x"], 'y': at0["y"], 'z': at0["z"]}
                    , 'end': {'x': at1["x"], 'y': at1["y"], 'z': at1["z"]}
                    , **kwargs})
    
def draw_sphere(view, at, kwargs={}):
    view.addSphere({'center': {'x': at['x'], 'y':at['y'], 'z':at['z']}, **kwargs })
