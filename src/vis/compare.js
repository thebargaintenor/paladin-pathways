/**
 * Ted Monchamp
 * MCSB913 - Applied Bioinformatics
 * Pathway Comparison Visualization
 *
 * --
 *
 * compare.js
 *   Visualization rendering script
 *
 * requires:
 *   index.html
 *   application.css
 *   manifest.json (and contents given in manifest)
 *   d3.min.js
 *   palette.js
 *
 * --
 *
 * This script leverages the D3.js framework to rebuild the pathway
 * diagram described in KEGG's KGML files, which have been converted
 * to JSON for ease of loading into browser for rendering.
 *
 * In order to keep from needing to regenerate this script for
 * different pathways, all file/data locations are placed in
 * manifest.json, in this same directory.  The manifest is loaded
 * at script startup and points this script toward all data sources.
 *
 * Created by tmonchamp on 4/18/16.
 */

// http://bl.ocks.org/dustinlarimer/5888271
var marker_paths = [
    { id: 0, name: 'circle', path: 'M 0, 0  m -5, 0  a 5,5 0 1,0 10,0  a 5,5 0 1,0 -10,0', viewbox: '-6 -6 12 12' },
    { id: 1, name: 'square', path: 'M 0,0 m -5,-5 L 5,-5 L 5,5 L -5,5 Z', viewbox: '-5 -5 10 10' },
    { id: 2, name: 'arrow', path: 'M 0,0 m -5,-5 L 5,0 L -5,5 Z', viewbox: '-5 -5 10 10' },
    { id: 2, name: 'stub', path: 'M 0,0 m -1,-5 L 1,-5 L 1,5 L -1,5 Z', viewbox: '-1 -5 2 10' }
];

var manifest;
var dataset;
var pathway_data;
var compound_hash = {};

var diagram = d3.select('#diagram');

var label_margin = 5;

var tooltip = d3.select('body').append('div')
    .attr('class', 'tooltip')
    .style('opacity', 0);

var defs = diagram.append('svg:defs');


// main items
d3.json('manifest.json', function(error, m) {
    if (error) {
        console.log(error);
    }
    else {
        manifest = m;
        console.log(manifest);

        d3.json(manifest['src'], function(error, d) {
            if (error) {
                console.log(error);
            }
            else {
                loadPathway(d);
            }
        });
    }
});

function generateMarkers(palette) {
    var arrow = marker_paths[2];
    defs.selectAll('marker')
        .data(palette)
        .enter()
        .append('svg:marker')
        .attr({
            'id' : function(c) { return 'marker_' + c; },
            'markerHeight': 10,
            'markerWidth': 10,
            'markerUnits': 'strokeWidth',
            'orient': 'auto',
            'refX': 0,
            'refY': 0,
            'viewBox': function() { return arrow.viewbox; }
        })
        .append('svg:path')
        .attr({
            'd': function() { return arrow.path; },
            'fill' : function(c) { return '#' + c; }
        });
}

function loadPathway(data) {
    pathway_data = data['pathway'];

    document.title += ' - ' + pathway_data['@title'];

    d3.select('#pathway_name').text(pathway_data['@title']);
    d3.select('#pathway_id').text('KEGG ID: ' + pathway_data['@number']);

    if (pathway_data['@description']) {
        d3.select('#pathway_description').text(pathway_data['@description']);
    }

    // process paths and draw items as needed
    var compounds = pathway_data['entry'];

    // build hash table to hold IDs and names
    compounds.forEach(function (c) {
        // now implementing ability to have a list of ECs as dependency
        compound_hash[c['@id']] = c['@name'].split(' ').map(function(ec) { return ec.substring(3); });
    });

    // add the group for paths first
    var paths = diagram.append('svg:g').attr({id : 'paths'});

    diagram.selectAll('enzyme')
        .data(compounds.filter(function(c) {return c['@type'] == 'enzyme';}))
        .enter()
        .append('rect')
        .attr({ // NOTE THAT KGML COORDINATES ARE ALL CENTER POINTS
            x: function (c) { return parseInt(c['graphics']['@x']) - parseInt(c['graphics']['@width']) / 2; },
            y: function (c) { return parseInt(c['graphics']['@y']) - parseInt(c['graphics']['@height']) / 2; },
            width: function (c) { return parseInt(c['graphics']['@width']); },
            height: function (c) { return parseInt(c['graphics']['@height']); },
            'class': 'enzyme',
            id: function(c) { return 'enzyme-' + c['@id']; }
        })
        .on('mouseover', function(d) {
            var contents = d['graphics']['@name'];
            if (manifest['dataset'] && pathway_data['enzymes']) {
                contents = buildTooltip(compound_hash[d['@id']]); // (d['graphics']['@name'])
            }

            if (contents) {
                tooltip
                    .style({
                        'opacity': 0.9,
                        'left': (d3.event.pageX) + 'px',
                        'top': (d3.event.pageY) + 'px'
                    })
                    .html(contents);
            }
        })
        .on('mouseout', function(d) {
            tooltip.style('opacity', 0);
        });

    // label enzymes in pathway
    diagram.selectAll('enzyme-label')
        .data(compounds.filter(function(c) {return c['@type'] == 'enzyme';}))
        .enter()
        .append('text')
        .attr({
            x: function (c) { return parseInt(c['graphics']['@x']); },
            y: function (c) {
                return parseInt(c['graphics']['@y']) + parseInt(c['graphics']['@height']) / 2 - label_margin;
            },
            'class': 'enzyme-label',
            'text-anchor': 'middle',
            id: function(c) { return 'enzyme-label-' + c['@id']; }
        })
        .text(function(c) { return c['graphics']['@name']; });
        //.text(function(c) { return c['@id']; });

    // label substrates/products on pathway
    diagram.selectAll('compound')
        .data(compounds.filter(function(c) {return c['@type'] == 'compound';}))
        .enter()
        .append('circle')
        .attr({
            cx: function (c) { return parseInt(c['graphics']['@x']); },
            cy: function (c) { return parseInt(c['graphics']['@y']); },
            r: function (c) { return parseInt(c['graphics']['@width']) / 2; },
            'class': 'compound',
            id: function(c) { return 'compound-' + c['@id']; }
        })
        .on('mouseover', function(d) {
            tooltip
                .style({
                    'opacity' : 0.9,
                    'left' : (d3.event.pageX) + 'px',
                    'top' : (d3.event.pageY) + 'px'
                })
                .html(pathway_data['compounds'][d['graphics']['@name']]); // + '<br/>' + d['@id']);
                //.html(d['graphics']['@name'] + '<br/>' +  pathway_data['compounds'][d['graphics']['@name']]);
        })
        .on('mouseout', function(d) {
            tooltip.style('opacity', 0);
        });

    diagram.selectAll('map')
        .data(compounds.filter(function(c) {
            return c['@type'] == 'map' && !c['graphics']['@name'].startsWith('TITLE:');
        }))
        .enter()
        .append('rect')
        .attr({ // NOTE THAT KGML COORDINATES ARE ALL CENTER POINTS
            x: function (c) { return parseInt(c['graphics']['@x']) - parseInt(c['graphics']['@width']) / 2; },
            y: function (c) { return parseInt(c['graphics']['@y']) - parseInt(c['graphics']['@height']) / 2; },
            rx: 10,
            ry: 10,
            width: function (c) { return parseInt(c['graphics']['@width']); },
            height: function (c) { return parseInt(c['graphics']['@height']); },
            'class': 'map',
            id: function(c) { return 'map-' + c['@id']; }
        });

    diagram.selectAll('map-label')
        .data(compounds.filter(function(c) {
            return c['@type'] == 'map' && !c['graphics']['@name'].startsWith('TITLE:');
        }))
        .enter()
        .append('text')
        .attr({
            x: function (c) { return parseInt(c['graphics']['@x']); },
            y: function (c) {
                return parseInt(c['graphics']['@y']) - parseInt(c['graphics']['@height']) / 2;
            },
            width: function (c) { return parseInt(c['graphics']['@width']); },
            'class': 'map-label',
            'text-anchor': 'middle',
            id: function(c) { return 'map-label-' + c['@id']; }
        })
        .text(function(c) { return c['graphics']['@name']; })
        .call(wrap);

    if(manifest['dataset']) {
        // fetch dataset file and draw links based on present enzymes
        d3.csv(manifest['dataset'], function (error, d) {
           if (error) {
               console.log(error);
           }
           else {
               drawPathways(d);
           }
        });
    }
    else {
        // draw full map with links
        // link behavior in KEGG has 4 options:
        // 1. vertical
        // 2. horizontal
        // 3. if sub.x < enz.x, change path x 
        paths.selectAll('default')
            .data(pathway_data['reaction'])
            .enter()
            .append('path')
            .attr({
                d: function (r) { return generateEnzymePath(r, 1, 1); },
                'stroke-width': 2,
                'stroke': 'white',
                'class': 'reaction default',
                'fill': 'transparent',
                id: function(r) { return 'reaction-' + r['@id']; }
            }); // find path attributes
    }

    // draw maplinks
    paths.selectAll('maplink')
        .data(pathway_data['relation'].filter(function (r) { return r['@type'] == 'maplink'; }))
        .enter()
        .append('path')
        .attr({
            'stroke-width': 2,
            'stroke': '#ccc',
            'class': 'maplink',
            'fill': 'transparent',
            'stroke-dasharray': '5 5',
            id: function (r) { return 'maplink-' + r['@entry1'] + '-' + r['@entry2']; },
            d: function (r) { return generatePath(r); }
        });
}

function buildTooltip(ec) {
    return dataset
        .map(function (set) { 
            // TODO: think about how we may want to represent compound entries
            return pathway_data['enzymes'][set['name']][ec[0]][0]; 
        })
        .filter(unique)
        .reduce(function (a, b) {
            if (!a) {
                return b;
            } else {
                return a + '<br/>' + b;
            }
        });
}

function drawPathways(d) {
    dataset = d;

    var paths = diagram.select('#paths');
    var legend = d3.select('#legend');

    var color_scheme = palette('qualitative', dataset.length, 0);
    generateMarkers(color_scheme);
    
    // mark gaps in the pathway
    paths.selectAll('default')
        .data(pathway_data['reaction'].filter(function (r) {
            var ecs = compound_hash[r['@id']]; // it's a list now!
            return dataset.every(function (set) {
                return !ecs.some(function(ec) { 
                    return !set.hasOwnProperty(ec) || set[ec] == 0; 
                });
            });
        }))
        //.data(pathway_data['reaction'])
        .enter()
        .append('path')
        .attr({
            d: function (r) { return generateEnzymePath(r, 1, 1); },
            'stroke-width': 1,
            'stroke': '#aaa',
            'class': 'reaction default',
            'fill': 'transparent',
            'stroke-dasharray': '2 4',
            id: function (r) { return 'reaction-' + r['@id']; }
        });

    // draw
    dataset.forEach(function (set, i) {
        var name = set['name'];
        paths.selectAll('.' + name)
            .data(pathway_data['reaction'].filter(function (r) {
                var ecs = compound_hash[r['@id']];
                return ecs.every(function(ec) {
                    return set.hasOwnProperty(ec) && set[ec] > 0;
                });
            }))
            .enter()
            .append('path')
            .attr({
                d: function (r) { return generateEnzymePath(r); },
                'stroke-width': 1,
                'stroke': function (r) { return '#' + color_scheme[i]; },
                'fill': 'transparent',
                'marker-end': function() { return 'url(#marker_' + color_scheme[i] + ')'; },
                'transform': function (r) { 
                    var offset = dataset.length / 2 - (i);
                    return 'translate(' + (-1 * offset) + ',' + offset + ')';
                },
                'class': function (r) { return 'reaction ' + name; },
                id: function(r) { return 'reaction-' + name + '-' + r['@id']; }
            });

        // add label to legend

        var label = legend.append('p');

        label.append('input')
            .attr({
                'type' : 'checkbox',
                'class' : 'series-toggle',
                'checked' : '1',
                'value' : name
            });

        label.append('div')
            .attr('class', 'legend-indicator')
            .style('background', '#' + color_scheme[i]);

        completeness = parseFloat(set['completeness']).toFixed(2);
        label.append('span').text(name + ' (' + completeness + ')');
    });

    d3.selectAll('.series-toggle').on('change', function () {
        var name = this.value;
        paths.selectAll('.' + name).style('opacity', this.checked ? 1 : 0);
    });
}

// I was hoping this could be more elegant, but this won't happen without ignoring
// all of the KEGG layout attributes.  Maybe later.
function generateEnzymePath(reaction, n, ofn) {
    var enzyme = d3.select('#enzyme-' + reaction['@id']);
    var substrate = d3.select('#compound-' + reaction['substrate']['@id']);
    var product = d3.select('#compound-' + reaction['product']['@id']);

    // ahahahahahahahaha I hate you javascript.
    var ex = parseFloat(enzyme.attr('x'));
    var ey = parseFloat(enzyme.attr('y'));
    var ew = parseFloat(enzyme.attr('width'));
    var eh = parseFloat(enzyme.attr('height'));
    var sx = parseFloat(substrate.attr('cx'));
    var sy = parseFloat(substrate.attr('cy'));
    var px = parseFloat(product.attr('cx'));
    var py = parseFloat(product.attr('cy'));

    var offset = 5;

    var dx = px - sx;
    var dirx = dx / Math.abs(dx);
    var dy = py - sy;
    var diry = dy / Math.abs(dy);

    var path = 'M' + sx + ',' + sy + ' ';
    if (sx == px) {
        // vertical
        path += 'V' + py;
    }
    else if (sy == py) {
        // horizontal
        path += 'H' + px;
    }
    else {
        // regrettably there's no string formatting in javascript.  rip.
        //console.log('Sy: ' + sy + ' Ey: ' + ey + ' Ey+Eh: ' + (ey + eh) + ' Total: ' + (sy >= ey && sy < (ey + eh)));

        var r = 15;

        if (inRange(sy, ey, eh) || inRange(px, ex, ew)) { // dx first
            //path += 'C' + px + ',' + sy + ' '  + px + ',' + sy + ' ' + px + ',' + py;
            path += 'h ' + (dx - dirx * r) + ' ';
            path += 'C ' + px + ',' + sy + ' ' + px + ',' + sy + ' ' + px + ',' + (sy + diry * r) + ' ';
            path += 'v ' + (dy - diry * r);
        }
        else if (inRange(py, ey, eh) || inRange(sx, ex, ew)) { // dy first
            // path += 'C' + sx + ',' + py + ' '  + sx + ',' + py + ' ' + px + ',' + py;
            path += 'v ' + (dy - diry * r) + ' ';
            path += 'C ' + sx + ',' + py + ' ' + sx + ',' + py + ' ' + (sx + dirx * r) + ',' + py + ' ';
            path += 'h ' + (dx - dirx * r);
        }
        else { // two bends: dy, dx, dy - by arbitrary designation
            var ecy = ey + eh / 2;
            var dcy = py - ecy;
            var dircy = dcy / Math.abs(dcy);

            path += 'v ' + (ecy - sy - diry * r) + ' ';
            path += 'C ' + sx + ',' + ecy + ' ' + sx + ',' + ecy + ' ' + (sx + dirx * r) + ',' + ecy + ' ';
            path += 'h ' + (dx - dirx * r * 2);
            path += 'C ' + px + ',' + ecy + ' ' + px + ',' + ecy + ' ' + px + ',' + (ecy + dircy * r) + ' ';
            path += 'v ' + (py - ecy - dircy * r) + ' ';
        }
    }

    return path;
}

function inRange(v, lbound, range) {
    return (v >= lbound && v < (lbound + range));
}

function generatePath(reaction) {
    var substrate = diagram.select('#compound-' + reaction['subtype']['@value']);
    var sx = parseFloat(substrate.attr('cx'));
    var sy = parseFloat(substrate.attr('cy'));

    // assume it's a maplink for now... can't find the compound to compound links yet
    var product = diagram.select('#map-' + reaction['@entry1']);
    if (product.empty()) {
        product = diagram.select('#map-' + reaction['@entry2']);
    }

    var px = parseFloat(product.attr('x'));
    var py = parseFloat(product.attr('y'));
    var pw = parseFloat(product.attr('width'));
    var ph = parseFloat(product.attr('height'));
    px += pw / 2;
    py += ph / 2;

    var r = 15;

    var dx = px - sx;
    var dirx = dx / Math.abs(dx);
    var dy = py - sy;
    var diry = dy / Math.abs(dy);

    var path = 'M' + sx + ',' + sy + ' ';
    if (sx == px) {
        // vertical
        path += 'V' + py;
    }
    else if (sy == py) {
        // horizontal
        path += 'H' + px;
    }
    else {
        // just do dx first
        if (px > sx) {
            //path += 'C' + px + ',' + sy + ' ' + px + ',' + sy + ' ' + px + ',' + py;
            path += 'h ' + (dx - dirx * r) + ' ';
            path += 'C ' + px + ',' + sy + ' ' + px + ',' + sy + ' ' + px + ',' + (sy + diry * r) + ' ';
            path += 'v ' + (dy - diry * r);
        }
        else {
            //path += 'C' + sx + ',' + py + ' '  + sx + ',' + py + ' ' + px + ',' + py;
            path += 'v ' + (dy - diry * r) + ' ';
            path += 'C ' + sx + ',' + py + ' ' + sx + ',' + py + ' ' + (sx + dirx * r) + ',' + py + ' ';
            path += 'h ' + (dx - dirx * r);
        }
    }

    return path;
}

// shamelessly pilfered from https://bl.ocks.org/mbostock/7555321
function wrap(text) {
    text.each(function() {
        var text = d3.select(this),
            words = text.text().split(/\s+/).reverse(),
            word,
            line = [],
            lineNumber = 1,
            lineHeight = 1.1, // ems
            width = text.attr('width'),
            x = text.attr('x'),
            y = text.attr('y'),
            tspan = text.text(null).append('tspan').attr({
                x: x,
                y: y,
                dy: lineHeight + 'em'
            });
        while (word = words.pop()) {
            line.push(word);
            tspan.text(line.join(' '));
            if (tspan.node().getComputedTextLength() > width) {
                line.pop();
                tspan.text(line.join(' '));
                line = [word];
                tspan = text.append('tspan').attr({
                    x: x,
                    y: y,
                    'dy': ++lineNumber * lineHeight + 'em'
                })
                .text(word);
            }
        }
    });
}

function unique(value, index, self) {
    return self.indexOf(value) === index;
}