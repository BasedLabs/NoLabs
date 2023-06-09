var localisation = function (probabilities) {
    window.onmousemove = function (e) {
        //console.log("mouse location:", e.clientX, e.clientY)
    };
    if (!document.getElementById('img'))
        throw Error('Declare img tag');
    const _createImage = (src) => {
        const image = new Image(500, 500);
        image.src = src;
        return image;
    };
    const img = document.getElementById('img');
    const obj = {
        hoverControlHighlight: {
            onmouseover: (el) => {
                const jqueryEl = $(el.target);
                jqueryEl.addClass('active');
                const localisationObjectKey = jqueryEl.data('localisationObjectKey');
                const highlightElement = obj.highlightData[localisationObjectKey];
                img.src = highlightElement.image.src;
            },
            onmouseleave: (el) => {
                $(el.target).removeClass('active');
                img.src = obj.highlightData.original.image.src;
            },
        },
        highlightData: {
            mithochondria: {
                image: _createImage('./static/images/mithochondria.png'),
                controlElementId: 'mithochondria-list-item',
                text: 'Mithochondria',
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
                controlElementId: 'nucleus-list-item',
                text: 'Nucleus',
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
                controlElementId: 'cytoplasm-list-item',
                text: 'Cytoplasm',
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
            other: {
                image: _createImage('./static/images/original.png'),
                controlElementId: 'other-proteins-item',
                text: 'Other proteins',
                polygon: [
                    [
                        [0, 0],
                        [500, 0],
                        [500, 500],
                        [0, 500]
                    ]
                ]
            },
            extracellular: {
                image: _createImage('./static/images/original.png'),
                controlElementId: 'extracellular-proteins-item',
                text: 'Extracellular proteins',
                polygon: [
                    [
                        [0, 0],
                        [500, 0],
                        [500, 500],
                        [0, 500]
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
            const highlightElements = [];
            for (const key in obj.highlightData) {
                const data = obj.highlightData[key];
                const elId = data.controlElementId;
                const el = $('#' + elId);
                el.append(data.text + ' ' + (probabilities[key].toFixed(2) * 100) + '%');
                el.hover(obj.hoverControlHighlight.onmouseover,
                    obj.hoverControlHighlight.onmouseleave);
                el.data('localisationObjectKey', key);
                highlightElements.push(el);
            }

            img.src = obj.highlightData.original.image.src;
            img.addEventListener('mousemove', event => {
                const rect = img.getBoundingClientRect();
                const x = event.clientX - rect.left;
                const y = event.clientY - rect.top;
                for (const key in obj.highlightData) {
                    const data = obj.highlightData[key];
                    if (obj.isInside([x, y], data.polygon)) {
                        for (const highlightElement of highlightElements) {
                            highlightElement.removeClass('active');
                        }
                        $('#' + data.controlElementId).addClass('active');
                        img.src = data.image.src;
                        return;
                    }
                }
            });
            img.addEventListener('mouseout', function () {
                img.src = obj.highlightData.original.image.src;
                for (const highlightElement of highlightElements) {
                    highlightElement.removeClass('active');
                }
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

