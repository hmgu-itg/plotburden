showhide_sp_code="""
    alert("hello")
    var d=source.data
    var alpha=d.alpha
    if(cb_obj.active==0){
    // hide sp
        for (var k = 0, len = alpha.length; k < alpha.length; k++){
        if(alpha[k]==0.7){alpha[k]=0}
    }
    }else{
    //show sp
    for (var k = 0, len = alpha.length; k < alpha.length; k++){
        if(alpha[k]==0){alpha[k]=0.7}
    }
    }
    d.alpha=alpha
    source.change.emit()
    """
changecolor_code="""
    d=source.data
    color=d.color
    if(cb_obj.active==0){
    // no specific color
        for (k = 0, len = color.length; k < color.length; k++){
        color[k]="#1F77B4"
    }
    }else if (cb_obj.active==1){
    //show sp
    for (k = 0, len = color.length; k < color.length; k++){
        color[k]=d.mafcolor[k]
    }
    }else{
        for (k = 0, len = color.length; k < color.length; k++){
        color[k]=d.weightcolor[k]
        }
    }
    d.color=color
    source.change.emit()
    """

changecohort_code="""
    console.log(names.data)
    console.log(names.data.toString())
    var co_names=names.data.split(',')
    var source=mainsource.data

    var f = cb_obj.value;
    for(k=0, len=co_names.length, k<co_names.length;k++){
        if(f==co_names[k]){
            source=eval(co_names[k])
            source.change.emit()
        }
    }
"""

hideburden_code="""
    d=source.data
    alpha=d.alpha
    if(cb_obj.active==0){
    // hide
        alpha[0]=1

    }else{
    //show
            alpha[0]=0
    }
    d.alpha=alpha
    source.change.emit()
    """

ld_hover_code="""
var lddata = lds.data;
var rdat = rawdat.data;
var c=["#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#ffffbf", "#fee08b", "#fdae61", "#f46d43"," #d53e4f"];
    console.log(cb_data.index)

var indices = cb_data.index.indices;
if(indices.length>0){
    for (var k = 0, len = indices.length; k < indices.length; k++){
        for (var i = 0; i < rdat['ps'].length; i++) {
            rdat['outalpha'][i]=0
            rdat['outcol'][i]="#3288bd"
        }
        burdenpos=rdat['ps'][indices[k]]
        for (var i = 0; i < rdat['ps'].length; i++) {
            if(rdat['alpha'][i]>0){
                if(i==indices[k]) {
                    rdat['outcol'][i]=c[8]
                    rdat['outalpha'][i]=1
                }
                else {
                rdat['outcol'][i]=c[0];
                    for (var j = 0; j < lddata['x1'].length; j++) {
                        if(lddata['x1'][j]==rdat['ps'][i] && lddata['x2'][j]==burdenpos){
                            rdat['outalpha'][i]=1
                            var r2=lddata['r2'][j];
                            switch(true){
                                case(r2>0.9):
                                rdat['outcol'][i]=c[8]
                                break;
                                case(r2>0.8):
                                rdat['outcol'][i]=c[7]
                                break;
                                case(r2>0.7):
                                rdat['outcol'][i]=c[6]
                                break;
                                case(r2>0.6):
                                rdat['outcol'][i]=c[5]
                                break;
                                case(r2>0.5):
                                rdat['outcol'][i]=c[4]
                                break;
                                case(r2>0.4):
                                rdat['outcol'][i]=c[3]
                                break;
                                case(r2>0.3):
                                rdat['outcol'][i]=c[2]
                                break;
                                case(r2>0.2):
                                rdat['outcol'][i]=c[1]
                                break;
                                default:
                                rdat['outcol'][i]=c[9]
                            }
                        }
                        else if(lddata['x1'][j]==burdenpos && lddata['x2'][j]==rdat['ps'][i]){
                            rdat['outalpha'][i]=1
                            var r2=lddata['r2'][j];
                            switch(true){
                                case(r2>0.9):
                                rdat['outcol'][i]=c[8]
                                break;
                                case(r2>0.8):
                                rdat['outcol'][i]=c[7]
                                break;
                                case(r2>0.7):
                                rdat['outcol'][i]=c[6]
                                break;
                                case(r2>0.6):
                                rdat['outcol'][i]=c[5]
                                break;
                                case(r2>0.5):
                                rdat['outcol'][i]=c[4]
                                break;
                                case(r2>0.4):
                                rdat['outcol'][i]=c[3]
                                break;
                                case(r2>0.3):
                                rdat['outcol'][i]=c[2]
                                break;
                                case(r2>0.2):
                                rdat['outcol'][i]=c[1]
                                break;
                                default:
                                rdat['outcol'][i]=c[9]
                            }
                        }

                    }

            }
            }
        }
    }
}
rawdat.change.emit()

"""

hover_test_code="""
    console.log(cb_data.index.indices)
"""

ldbz_hover_code="""
var sgl = signalling.data;
if (sgl['way'][0] == 0) {
    // the signalling dataframe informs us that the first item in the radiolist is selected (highlight points)
    var lddata = lds.data;
    var rdat = rawdat.data;
    var c = ["#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#ffffbf", "#fee08b", "#fdae61", "#f46d43", " #d53e4f"];
    console.log(cb_data.index)
    var indices = cb_data.index.indices;
    if (indices.length > 0) {
        // this bit of code was higher, but was causing a reinit at every mouse move.
        // Now the LD display sticks, and reinit is done when a point is hovered.
        for (var i = 0; i < rdat['ps'].length; i++) {
                    rdat['outalpha'][i] = 0
                    rdat['outcol'][i] = "#3288bd"
                }

        for (var k = 0, len = indices.length; k < indices.length; k++) {
            var burdenpos = rdat['ps'][indices[k]]
            for (var i = 0; i < rdat['ps'].length; i++) {
                if (rdat['alpha'][i] > 0) {
                    if (i == indices[k]) {
                        rdat['outcol'][i] = c[8]
                        rdat['outalpha'][i] = 1
                    } else {
                        rdat['outcol'][i] = c[0];
                        for (j = 0; j < lddata['x1'].length; j++) {
                            if (lddata['x1'][j] == rdat['ps'][i] && lddata['x2'][j] == burdenpos) {
                                rdat['outalpha'][i] = 1
                                var r2 = lddata['r2'][j];
                                switch (true) {
                                    case (r2 > 0.9):
                                        rdat['outcol'][i] = c[8]
                                        break;
                                    case (r2 > 0.8):
                                        rdat['outcol'][i] = c[7]
                                        break;
                                    case (r2 > 0.7):
                                        rdat['outcol'][i] = c[6]
                                        break;
                                    case (r2 > 0.6):
                                        rdat['outcol'][i] = c[5]
                                        break;
                                    case (r2 > 0.5):
                                        rdat['outcol'][i] = c[4]
                                        break;
                                    case (r2 > 0.4):
                                        rdat['outcol'][i] = c[3]
                                        break;
                                    case (r2 > 0.3):
                                        rdat['outcol'][i] = c[2]
                                        break;
                                    case (r2 > 0.2):
                                        rdat['outcol'][i] = c[1]
                                        break;
                                    default:
                                        rdat['outcol'][i] = c[9]
                                }
                            } else if (lddata['x1'][j] == burdenpos && lddata['x2'][j] == rdat['ps'][i]) {
                                rdat['outalpha'][i] = 1
                                r2 = lddata['r2'][j];
                                switch (true) {
                                    case (r2 > 0.9):
                                        rdat['outcol'][i] = c[8]
                                        break;
                                    case (r2 > 0.8):
                                        rdat['outcol'][i] = c[7]
                                        break;
                                    case (r2 > 0.7):
                                        rdat['outcol'][i] = c[6]
                                        break;
                                    case (r2 > 0.6):
                                        rdat['outcol'][i] = c[5]
                                        break;
                                    case (r2 > 0.5):
                                        rdat['outcol'][i] = c[4]
                                        break;
                                    case (r2 > 0.4):
                                        rdat['outcol'][i] = c[3]
                                        break;
                                    case (r2 > 0.3):
                                        rdat['outcol'][i] = c[2]
                                        break;
                                    case (r2 > 0.2):
                                        rdat['outcol'][i] = c[1]
                                        break;
                                    default:
                                        rdat['outcol'][i] = c[9]
                                }
                            }

                        }

                    }
                }
            }
        }
    }
    rawdat.change.emit()

} else if (sgl['way'][0] == 1) {
    // the signalling dataframe indicates that the user wants a fountain display
    var lddata = lds.data;
    var rdat = rawdat.data;
    var bzdata = bezier.data;
    var c = ["#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#ffffbf", "#fee08b", "#fdae61", "#f46d43", " #d53e4f"];
        console.log(cb_data.index)

    var indices = cb_data.index.indices;
    if (indices.length > 0) {
        // Idem here, for speed, we reinitialise only if another proper point is hovered.
        bzdata['x0'] = [];
        bzdata['x1'] = [];
        bzdata['y0'] = [];
        bzdata['y1'] = [];
        bzdata['cx0'] = [];
        bzdata['cx1'] = [];
        bzdata['cy0'] = [];
        bzdata['cy1'] = [];
        bzdata['col'] = [];
        for (var k = 0, len = indices.length; k < indices.length; k++) {
            var burdenpos = rdat['ps'][indices[k]]
            var spval = rdat['logsp'][indices[k]]
            for (var i = 0; i < rdat['ps'].length; i++) {
                if (rdat['alpha'][i] > 0) {
                    for (var j = 0; j < lddata['x1'].length; j++) {
                        if (lddata['x1'][j] == rdat['ps'][i] && lddata['x2'][j] == burdenpos) {
                            var r2 = lddata['r2'][j];
                            if (r2 > 0.1) {
                                bzdata['x0'][i] = burdenpos;
                                bzdata['cx0'][i] = burdenpos;
                                bzdata['x1'][i] = rdat['ps'][i];
                                bzdata['cx1'][i] = rdat['ps'][i];
                                bzdata['y1'][i] = rdat['logsp'][i];
                                bzdata['cy1'][i] = rdat['logsp'][i] + 1;
                                bzdata['y0'][i] = spval;
                                bzdata['cy0'][i] = spval + 1;
                            }
                            switch (true) {
                                case (r2 > 0.9):
                                    bzdata['col'][i] = c[8]
                                    break;
                                case (r2 > 0.8):
                                    bzdata['col'][i] = c[7]
                                    break;
                                case (r2 > 0.7):
                                    bzdata['col'][i] = c[6]
                                    break;
                                case (r2 > 0.6):
                                    bzdata['col'][i] = c[5]
                                    break;
                                case (r2 > 0.5):
                                    bzdata['col'][i] = c[4]
                                    break;
                                case (r2 > 0.4):
                                    bzdata['col'][i] = c[3]
                                    break;
                                case (r2 > 0.3):
                                    bzdata['col'][i] = c[2]
                                    break;
                                case (r2 > 0.2):
                                    bzdata['col'][i] = c[1]
                                    break;
                                default:
                                    bzdata['col'][i] = c[9]
                            }
                        } else if (lddata['x1'][j] == burdenpos && lddata['x2'][j] == rdat['ps'][i]) {
                            var r2 = lddata['r2'][j];
                            if (r2 > 0.1) {
                                bzdata['x0'][i] = burdenpos;
                                bzdata['cx0'][i] = burdenpos;
                                bzdata['x1'][i] = rdat['ps'][i];
                                bzdata['cx1'][i] = rdat['ps'][i];
                                bzdata['y1'][i] = rdat['logsp'][i];
                                bzdata['cy1'][i] = rdat['logsp'][i] + 1;
                                bzdata['y0'][i] = spval;
                                bzdata['cy0'][i] = spval + 1;
                            }
                            switch (true) {
                                case (r2 > 0.9):
                                    bzdata['col'][i] = c[8]
                                    break;
                                case (r2 > 0.8):
                                    bzdata['col'][i] = c[7]
                                    break;
                                case (r2 > 0.7):
                                    bzdata['col'][i] = c[6]
                                    break;
                                case (r2 > 0.6):
                                    bzdata['col'][i] = c[5]
                                    break;
                                case (r2 > 0.5):
                                    bzdata['col'][i] = c[4]
                                    break;
                                case (r2 > 0.4):
                                    bzdata['col'][i] = c[3]
                                    break;
                                case (r2 > 0.3):
                                    bzdata['col'][i] = c[2]
                                    break;
                                case (r2 > 0.2):
                                    bzdata['col'][i] = c[1]
                                    break;
                                default:
                                    bzdata['col'][i] = c[9]
                            }
                        }

                    }


                }
            }
        }
    }
    bezier.change.emit()
}
"""

changehover_code="""
var sgl=signalling.data;
if(cb_obj.active==0){
    // used to be fountain, is now highlight. Destroy fountain.
    var bzdata = bezier.data;
    bzdata['x0'] = [];
    bzdata['x1'] = [];
    bzdata['y0'] = [];
    bzdata['y1'] = [];
    bzdata['cx0'] = [];
    bzdata['cx1'] = [];
    bzdata['cy0'] = [];
    bzdata['cy1'] = [];
    bzdata['col'] = [];
    bezier.change.emit();
}else if (cb_obj.active==1){
    // used to be highlight, is now fountain. Destroy highlight.
        var rdat = rawdat.data;
    for (i = 0; i < rdat['ps'].length; i++) {
            rdat['outalpha'][i] = 0
            rdat['outcol'][i] = "#3288bd"
    }
    rawdat.change.emit()
}
sgl['way'][0]=cb_obj.active;
signalling.change.emit()
"""

displayhits_code="""
    d=source.data
    alpha=d['alpha']
    if(cb_obj.active==0){
    // hide signals
        for (k = 0, len = alpha.length; k < alpha.length; k++){
        alpha[k]=0
    }
    }else if (cb_obj.active==1){
    //show signals
        for (k = 0, len = alpha.length; k < alpha.length; k++){
            if(d['pheno']!="none"){
                alpha[k]=1
            }
        }
    }
    d['alpha']=alpha
    source.change.emit()
    """
