var localisation = function (probabilities) {
    if (!document.getElementById('img'))
        throw Error('Declare img tag');
    const _createImage = (src) => {
        const image = new Image(500, 500);
        image.src = src;
        return image;
    };
    const obj = {
        highlightData: {
            mithochondria: {
                image: _createImage('./static/images/mithochondria.png'),
                polygon: [
                    [
                        [150, 115],
                        [182, 112],
                        [217, 139],
                        [210, 155],
                        [180, 143],
                        [156, 149]
                    ],
                    [
                        [213, 352],
                        [234, 345],
                        [250, 353],
                        [242, 385],
                        [200, 390],
                        [182, 377],
                        [178, 353],
                        [192, 345],
                    ]
                ]
            },
            nucleus: {
                image: _createImage('./static/images/nucleus.png'),
                polygon: [
                    [[176, 178],
                        [202, 199],
                        [222, 229],
                        [216, 257],
                        [170, 271],
                        [134, 255],
                        [137, 201],
                        [169, 177]]
                ]
            },
            cytoplasm: {
                image: _createImage('./static/images/cytoplasm.png'),
                polygon: [
                    [
                        [95, 160],
                        [190, 99],
                        [238, 167],
                        [341, 203],
                        [400, 313],
                        [332, 385],
                        [160, 381],
                        [88, 315],
                        [82, 201],
                        [130, 107]
                    ]
                ]
            },
            original: {
                image: _createImage('./static/images/original.png'),
                polygon: [
                    [
                        [0, 0],
                        [500, 0],
                        [500, 500],
                        [0, 500]
                    ]
                ]
            },
        },
        isInside: function inside(point, vsArray) {
            for (const vs of vsArray) {
                // ray-casting algorithm based on
                // https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html

                const x = point[0], y = point[1];

                let inside = false;
                for (let i = 0, j = vs.length - 1; i < vs.length; j = i++) {
                    let xi = vs[i][0], yi = vs[i][1];
                    let xj = vs[j][0], yj = vs[j][1];

                    let intersect = ((yi > y) != (yj > y))
                        && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
                    if (intersect) inside = !inside;
                }

                if (inside)
                    return true;
            }
            return false;
        },
        render: function () {
            const img = document.getElementById('img');
            img.src = obj.highlightData.original.image.src;
            img.addEventListener('mousemove', event => {
                const rect = img.getBoundingClientRect();
                const x = event.clientX - rect.left;
                const y = event.clientY - rect.top;
                for (const key in obj.highlightData) {
                    const data = obj.highlightData[key];
                    if (obj.isInside([x, y], data.polygon)) {
                        img.src = data.image.src;
                        return;
                    }
                }
            });
            img.addEventListener('mouseout', function () {
                img.src = obj.highlightData.original.image.src;
            });
        }
    }
    return obj;
}
var folding = function (proteinFileContent) {
    const stage = new NGL.Stage("viewport");
    stage.setParameters({backgroundColor: 'white'})
    const obj = {
        views: {
            default: 'default',
            cartoon: 'cartoon',
            backbone: 'backbone',
            ballsAndSticks: 'ball+stick',
            contact: 'contact',
            helixorient: 'helixorient',
            hyperball: 'hyperball',
            licorice: 'licorice',
            ribbon: 'ribbon',
            rope: 'rope',
            surface: 'surface',
            spacefill: 'spacefill',
            unitcell: 'unitcell'
        },
        reload: () => {
            stage.removeAllComponents();
            const proteinFileContentBlob = new Blob([proteinFileContent], {type: 'text/plain'});
            const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', {type: 'text/plain'});
            stage.loadFile(proteinFile, {defaultRepresentation: true}).then((component) => {
                obj.component = component;
            });
        },
        setView: (viewName) => {
            if (viewName === obj.views.default) {
                obj.reload();
                return;
            }
            obj.component.removeAllRepresentations();
            obj.component.addRepresentation(viewName);
        },
        render: () => {
            document.addEventListener("DOMContentLoaded", () => {
                obj.reload();
            });
        }
    }
    return obj;
}

